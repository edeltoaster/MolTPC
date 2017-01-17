/*
MolTPC framework
Copyright (C) 2011-2013  Thorsten Will

This program is free software; you can redistribute it and/or modify it under the terms of 
the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, 
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; 
if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110, USA
*/

package base;
import java.util.Collection;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentSkipListSet;

/**
 * Conformational search class
 * uses metropolis-criterion for traversal
 * works in parallel
**/
public class ConfSearch { 
	private int temperature = 295;
	protected Random rnd = new Random(System.currentTimeMillis());
	private static final double R = 8.314472;
	protected final ConfSearchTree tree_root;
	protected final double step;
	protected final int no_steps;
	protected final ConcurrentSkipListSet<String> fast_ring_closures = new ConcurrentSkipListSet<String>(); // thread-safe
	protected final ConcurrentSkipListSet<String> rings_smiles = new ConcurrentSkipListSet<String>(); // thread-safe
	protected final ConcurrentLinkedQueue<Molecule> new_found = new ConcurrentLinkedQueue<Molecule>(); // thread-safe
	protected final double E; //dielectric constant
	protected final Runtime rt = Runtime.getRuntime(); // for memory-management
	protected int n_DEE = 0;
	protected final boolean use_DEE;
	
	public ConfSearch(Molecule molecule, int no_steps, double E,int temperature,boolean dee) {
		this(molecule,no_steps,E,dee);
		this.temperature = temperature;
	}
	
	public ConfSearch(List<Molecule> molecules, int no_steps, double E,int temperature,boolean dee) {
		this(molecules,no_steps,E,dee);
		this.temperature = temperature;
	}
	
	public ConfSearch(Molecule molecule,int no_steps,double E, boolean dee) {
		this.step = 2*Math.PI/no_steps;
		this.tree_root = new ConfSearchTree(null,null);
		this.no_steps = no_steps;
		this.E = E;
		this.use_DEE = dee;
		List<Molecule> molecules = new LinkedList<Molecule>();
		molecules.add(molecule);
		build(molecules);
	}
	
	public ConfSearch(List<Molecule> molecules,int no_steps,double E,boolean dee) {
		this.step = 2*Math.PI/no_steps;
		this.tree_root = new ConfSearchTree(null,null);
		this.no_steps = no_steps;
		this.E = E;
		this.use_DEE = dee;
		build(molecules);
	}
	
	private void build(List<Molecule> molecules) {
		Deque<Molecule> stack = new LinkedList<Molecule>(molecules);
		
		RotationBuilder[] threads = null;
		
		// workin'
		while (stack.size()>0) {
		
			threads = new RotationBuilder[Math.min(stack.size(),BabelInterface.threads)];
		
			// do-loop for max utilization
			do {
				// start jobs
				for (int i=0;i<threads.length;i++) {
					threads[i] = new RotationBuilder(tree_root,stack.pop());
					threads[i].start();
				}
				
				try { // wait for them getting done
					
					for (int i=0;i<threads.length;i++)
						threads[i].join();
					
				} catch (Exception e) {
					System.out.println("Problem in ConfSearch:RotationBuilder running ...");
					e.printStackTrace();
				}
				
				// checking free memory after every run, if desired ...
				if (use_DEE)
					CheckMemory();
				
			} while (stack.size()>=BabelInterface.threads); // less refilling more working
			
			// adding molecules with ring closure
			if (new_found.size() > 0) {
				stack.addAll(BabelInterface.batch_sd_minimize(new_found,100,(float)E));
				new_found.clear();
			}
	
		}
		
		// done workin'
		//tree_root.DEE(); //probably better with size-check?
		tree_root.recomputeBorders(false);// computes root interval
		tree_root.resetVisits();
		cleanup();
		
	}
	
	public void minimizeTree() {
		int n_old = tree_root.size();
		int n_new = n_old+1;
		
		while (n_old != n_new) {
			n_old = n_new;
			tree_root.DEE();
			n_new = tree_root.size();
		}
	}
	
	/**
	 * 
	 * worker class
	 *
	 */
	
	class RotationBuilder extends Thread {
		
		private ConfSearchTree parent;
		private Molecule molecule;
		
		RotationBuilder(ConfSearchTree parent,Molecule molecule) {
			this.parent = parent;
			this.molecule = molecule;
		}
		
		@Override
		public void run() {
			Map<Bond,List<Atom>> rot_map = molecule.getRotatableBondsMap();
			Map<Bond,Set<Atom>> disj_atom_map = Molecule.getDisjunctRotAtomsMap(rot_map);
			
			Bond[] rot_bonds = new Bond[rot_map.keySet().size()];
			
			// if there is nothing to rotate -> done
			if (rot_bonds.length==0){
				new TerminalNode(parent, molecule,-1,null,null,null,0);
				return;
			}
			
			int i = 0;
			
			// bonds sorted by atoms to rotate size
			for (Bond b:rot_map.keySet())
				rot_bonds[i++] = b;
			
			Bond max,current;
			int max_ind;
			for (i=0;i<rot_bonds.length-1;i++) {
				max = rot_bonds[i];
				max_ind = i;
				for (int j=i+1;j<rot_bonds.length;j++) {
					current = rot_bonds[j];
					if (rot_map.get(current).size()>rot_map.get(max).size()) {
						max = current;
						max_ind = j;
					}
				}
				rot_bonds[max_ind] = rot_bonds[i];
				rot_bonds[i] = max;
			}
			
			// bonds now sorted and ready to use :-P
			double[] fragment_score = new double[rot_bonds.length];
			for (i=0;i<rot_bonds.length;i++)
				fragment_score[i] = getScore(E,molecule,disj_atom_map.get(rot_bonds[i]));
			
			new RotationNode(parent,molecule,rot_bonds,rot_map,disj_atom_map,0,-1,molecule.getFixAtoms(disj_atom_map),fragment_score);
			return;
		}
		
	}
	
	private void generateMolecule(Molecule molecule,Atom O1,Atom O2) {
		
		// giving every guy a name
		Atom H1 = O1.getAdjHs().get(0);
		
		// generate some fingerprint
		int first = O1.getNo();
		int sec = O2.getNo();
		if (O1.getNo()>O2.getNo()) {
			first = O2.getNo();
			sec = O1.getNo();
		}
				
		String id = molecule.getInstance()+" "+first+" "+sec;
		
		// if already there from earlier conformation -> throw away
		if (fast_ring_closures.contains(id))
			return;
		
		fast_ring_closures.add(id);
		
		// ring size constraint
		for (Atom x:O1.getAdjFrame(0,true).get(0).getAdjAtoms())
			for (Atom y:O2.getAdjFrame(0,true).get(0).getAdjAtoms()) {
				if (x==y)
					return;
			}
		
		Molecule temp_mol = new Molecule(molecule,O1,H1,O2,false);
		String smiles = BabelInterface.getSMILESfast(temp_mol);
		
		// if already there -> throw away
		if (!rings_smiles.contains(smiles)) {
			rings_smiles.add(smiles);
			new_found.add(temp_mol);
		}
		
		// same for stereoform
		temp_mol = new Molecule(molecule,O1,H1,O2,true); // variant
		smiles = BabelInterface.getSMILESfast(temp_mol);
		
		// needs to be checked separately
		if (!rings_smiles.contains(smiles)) {
			rings_smiles.add(smiles);
			new_found.add(temp_mol);
		}
	
	}
	
	
	protected double getScore(double E_0,Molecule mol,Collection<Atom> atoms_to_score,Collection<Atom> scored_atoms) {
		double score = 0;
		double r,t1,t2,diff;
		
		for (Atom a:atoms_to_score) {
			for (Atom b:mol.score_map.get(a)) {
				if (scored_atoms.contains(b)) {
					
					r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
							+(a.getY()-b.getY())*(a.getY()-b.getY())
							+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
					
					// clash and RC detection
					if (r<3.45f) { // 3.45 because of radii from O in OH for OH detect.
						diff = a.getVdW()+b.getVdW() - r;
						// clash
						if (diff>1.5f) 
							return Double.NaN;
						else if (diff>0) {
							if (a.isOH() && b.isKeto())
								generateMolecule(mol,a,b); // generate molecule with ring-closure
							else if (a.isKeto() && b.isOH())
								generateMolecule(mol,b,a); // generate molecule with ring-closure
						}
					}
					
					// VdW
					t1 = (a.getVdW()+ b.getVdW())/r;
					double term6 = t1 * t1 * t1;
					term6 = term6 * term6;
					double term12 = term6 * term6;
					
					t2 = Math.sqrt(a.getVdWEpsilon()*b.getVdWEpsilon());
					
					score += 0.5* t2 * (term12 - 2*term6);
					
					
					// coulomb in kcal / mol
					score += 331.92552431411053 * a.getCharge()*b.getCharge() / ( r*E_0);	
				}
			}
			
			// for 1-4
			for (Atom b:mol.score_map_4.get(a)) {
				if (scored_atoms.contains(b)) {
					
					r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
							+(a.getY()-b.getY())*(a.getY()-b.getY())
							+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
					
					// clash and RC detection
					if (r<3.45f) { // 3.45 because of radii from O in OH for OH detect.
						diff = a.getVdW()+b.getVdW() - r;
						// clash
						if (diff>1.5f) 
							return Double.NaN;
						else if (diff>0) {
							if (a.isOH() && b.isKeto())
								generateMolecule(mol,a,b); // generate molecule with ring-closure
							else if (a.isKeto() && b.isOH())
								generateMolecule(mol,b,a); // generate molecule with ring-closure
						}
					}
					
					// VdW
					t1 = (a.getVdW()+ b.getVdW())/r;
					double term6 = t1 * t1 * t1;
					term6 = term6 * term6;
					double term12 = term6 * term6;
					
					t2 = Math.sqrt(a.getVdWEpsilon()*b.getVdWEpsilon());
					
					score += 0.5* t2 * (term12 - 2*term6);
					
					
					// coulomb in kcal / mol
					score += 0.5 * 331.92552431411053 * a.getCharge()*b.getCharge() / ( r*E_0);	
				}
			}
		}
			
		return score;
	}
	
	// done on first placed fragment
	protected double getScore(double E_0,Molecule mol,Collection<Atom> atoms_to_score) {
		double score = 0;
		double r,t1,t2;
		
		if (mol.score_map==null)
			mol.buildFastScoringMap();
		
		for (Atom a:atoms_to_score) {
			
			for (Atom b:mol.score_map.get(a)) {
				if (atoms_to_score.contains(b)) {
					
					r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
							+(a.getY()-b.getY())*(a.getY()-b.getY())
							+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
					
					// VdW
					t1 = (a.getVdW()+ b.getVdW())/r;
					double term6 = t1 * t1 * t1;
					term6 = term6 * term6;
					double term12 = term6 * term6;
					
					t2 = Math.sqrt(a.getVdWEpsilon()*b.getVdWEpsilon());
					
					score += t2 * (term12 - 2*term6);
					
					
					// coulomb in kcal / mol
					score += 331.92552431411053 * a.getCharge()*b.getCharge() / ( r*E_0);	
				}
			}
			
			// for 1-4
			for (Atom b:mol.score_map_4.get(a)) {
				if (atoms_to_score.contains(b)) {
					
					r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
							+(a.getY()-b.getY())*(a.getY()-b.getY())
							+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
					
					// VdW
					t1 = (a.getVdW()+ b.getVdW())/r;
					double term6 = t1 * t1 * t1;
					term6 = term6 * term6;
					double term12 = term6 * term6;
					
					t2 = Math.sqrt(a.getVdWEpsilon()*b.getVdWEpsilon());
					
					score += 0.5 * t2 * (term12 - 2*term6);
					
					
					// coulomb in kcal / mol
					score += 0.5 * 331.92552431411053 * a.getCharge()*b.getCharge() / ( r*E_0);	
				}
			}
		}
			
		return score;
	}
	
	public ConfSearchTree getSearchTree() {
		return tree_root;
	}
	
	public int getNoFoundRingClosures() {
		return rings_smiles.size();
	}
	
	public int getNumberOfDEECalls() {
		return n_DEE;
	}
	
	public List<Molecule> getNBest(int n) {
		List<Molecule> best = new LinkedList<Molecule>();
		int current = 0;
		
		if (n>tree_root.size())
			n = tree_root.size();
		
		Molecule temp;
		while(current < n) {
			temp = getBest();
			if (temp==null)
				return best;
			best.add(temp);
			current++;
		}
		
		return best;
	}
	
	private Molecule getBest() {
		ConfSearchTree current = tree_root;
		double min = tree_root.min_score;
		ConfSearchTree best = null;
		ConfSearchTree alternative = null;
		
		for (ConfSearchTree t:current.children)
			if (t.min_score==min) {
				best = t;
				break;
			}
		
		// stochastics incoming !!
		int n = 0;
		int use = rnd.nextInt(current.children.size());
		for (ConfSearchTree t:current.children) {
			if (n == use) {
				alternative = t;
				break;
			}
			n++;
		}
		
		double dE = alternative.min_score-best.min_score;
		if (MetropolisCrit(dE))
			current = alternative;
		else
			current = best;
		// end stochastics, current set
		
		// first RotationNode or TerminalNode
		if (current.molecule==null)
			return null;
		
		Molecule molecule = new Molecule(current.molecule);
		Map<Bond,List<Atom>> rot_map = molecule.getRotatableBondsMap();
		
		while (!(current instanceof TerminalNode)) {
			min = current.min_score;
			RotationNode rot = (RotationNode) current;
			for (ConfSearchTree t:current.children)
				if (t.min_score==min) {
					best = t;
					break;
				}
			
			// stochastics incoming !! (the second)
			n = 0;
			use = rnd.nextInt(current.children.size());
			for (ConfSearchTree t:current.children) {
				if (n == use) {
					alternative = t;
					break;
				}
				n++;
			}
			
			dE = alternative.min_score-best.min_score;
			if (MetropolisCrit(dE))
				current = alternative;
			else
				current = best;
			// end stochastics, current set
			
			// rotate
			if (current instanceof TerminalNode)
				rot.rotate(((TerminalNode)current).step_before, rot_map);
			else
				rot.rotate(((RotationNode)current).step_before, rot_map);
		}
		
		// print scoring of current
		//System.out.println(current.min_score);
		
		// cut current
		current.prune();
		
		// update parent(s)
		current.getParent().recomputeBorders(true);
		
		return molecule;
	}
	
	public int getTemperature() {
		return temperature;
	}
	
	protected boolean MetropolisCrit(double dE) {
		double boltzmann_prob = Math.exp(-dE/(temperature * R /4184.0));
		return (boltzmann_prob >= rnd.nextDouble()); // take it nontheless :-)
	}
	
	private void cleanup() {
		BabelInterface.minimizeClean();
		BabelInterface.smiClean();
	}
	
	/**
	 * frees memory using an assumptionfree analogon of the DEE-criterion if necessary
	 */
	private void CheckMemory() {
		// when to start DEE
		 if(rt.freeMemory()/rt.totalMemory()<=0.5) {
			 // start DEE
			 n_DEE++;
			 tree_root.DEE();
			 // clean memory
			 rt.gc();
		 }
	}
	
	/**
	 * Enforces n calls of the DEE-function
	 * @param n number of DEE-calls
	 */
	public void enforceDEE(int n) {
		for (int i=0;i<n;i++)
			tree_root.DEE();
	}
	
	
	/** specifies interface for traversal and pruning etcpp **/
	public class ConfSearchTree {
		protected final ConfSearchTree parent;
		protected final Molecule molecule;
		protected final ConcurrentLinkedQueue<ConfSearchTree> children = new ConcurrentLinkedQueue<ConfSearchTree>(); // thread-safe
		protected double min_score = Double.POSITIVE_INFINITY;
		protected double max_score = Double.NEGATIVE_INFINITY;
		protected double part_score = Double.POSITIVE_INFINITY;
		protected int visits = 0; // used for dynamic pruning and possible broader traversal
		
		public ConfSearchTree(ConfSearchTree parent,Molecule molecule) {
			this.parent = parent;
			this.molecule = molecule;
			
			if (parent!=null)
				parent.children.add(this);
		}
		
		public int size() {
			int sum = 0;
			
			for (ConfSearchTree t:children)
				sum += t.size();
			
			return sum;
		}
		
		public double getMin_score() {
			return min_score;
		}
		
		public void setMin_score(double min_score) {
			this.min_score = min_score;
		}
		
		public double getMax_score() {
			return max_score;
		}
		
		public void setMax_score(double max_score) {
			this.max_score = max_score;
			
		}
		
		public ConfSearchTree getParent() {
			return parent;
		}
		
		public ConcurrentLinkedQueue<ConfSearchTree> getChildren() {
			return children;
		}
		
		public void prune() {
			parent.children.remove(this);
		}
		
		protected void recomputeBorders(boolean recursive) {
			if (recursive) {
				min_score = Double.POSITIVE_INFINITY;
				max_score = Double.NEGATIVE_INFINITY;
			}
			
			for (ConfSearchTree t:children) {
				if (t.min_score < min_score)
					min_score = t.min_score;
				
				if (t.max_score > max_score)
					max_score = t.max_score;
			}
			
			if (min_score == Double.POSITIVE_INFINITY && parent!=null) 
				prune();
			
			if (recursive && parent!=null)
				parent.recomputeBorders(true);
		}
		
		public void visited(int cut) {
			visits++;
			if (visits>=cut) {
				prune();
				parent.recomputeBorders(true);
			}
		}
		
		public int getVisits() {
			return visits;
		}
		
		public void resetVisits() {
			visits = 0;
			for (ConfSearchTree tr:children)
				tr.resetVisits();
		}
		
		public void DEE() {
			visits++;
			double biggest_min = Double.NEGATIVE_INFINITY;
			double smallest_max = Double.POSITIVE_INFINITY;
			int min_visits = Integer.MAX_VALUE;
			
			// finding the least DEE-visited nodes
			for (ConfSearchTree t:children) {
				if (t.getVisits()<min_visits)
					min_visits = t.getVisits();
			}
			
			for (ConfSearchTree t:children) {
				if (t.min_score>biggest_min)
					biggest_min = t.min_score;
				if (t.max_score<smallest_max)
					smallest_max = t.max_score;
			}
			
			// DEE crit
			if (biggest_min>smallest_max) {
				
				for (ConfSearchTree t:children) {
					
					if (t.min_score>smallest_max && t.min_score!=t.max_score) // DEE crit only on intervals
						t.prune();
				}
				
				recomputeBorders(true);
			}
			
			// to avoid DEEing some subtrees more often
			for (ConfSearchTree t:children) {
				if (t.getVisits() == min_visits)
					t.DEE();
			}
		}

	}
	
	
	class RotationNode extends ConfSearchTree {
		protected Quaternion rot_quat;
		protected Quaternion rot_quat_conj;
		protected float[] corr;
		protected int step_before;
		protected List<Atom> just_used;
		protected final Bond[] bonds;
		protected final int bond_index;
		
		public RotationNode(ConfSearchTree parent,Molecule molecule,Bond[] bonds,
				Map<Bond,List<Atom>> rot_map,Map<Bond,Set<Atom>> disj,
				int bond_index,int step_before,List<Atom> just_tested,double[] fragment_score) {
			super(parent,molecule);
			this.step_before = step_before;
			this.bonds = bonds;
			this.bond_index = bond_index;
			
			double temp;
			if (bond_index!=0) {
				temp = getScore(E,molecule,disj.get(bonds[bond_index-1]),just_tested);
				if (Double.isNaN(temp)) {
					prune();
					return;
				}
				part_score = parent.part_score + temp +fragment_score[bond_index-1];
			} else {
				part_score = getScore(E,molecule,just_tested);
			}
				
			int n = 0;
			
			/** building of the quaternions **/
			
			Atom c1 = bonds[bond_index].getA1();
			Atom c2 = bonds[bond_index].getA2();
			
			
			float x = c2.getX()-c1.getX();
			float y = c2.getY()-c1.getY();
			float z = c2.getZ()-c1.getZ();
			
			// for unit-quaternion
			float factor = (float) (1/Math.sqrt(x*x+y*y+z*z));
			
			float[] rot_vector = new float[]{factor*x,
					factor*y,
					factor*z};
			
			// the quaternion rotates with "step"-width in every iteration
			float sin_coeff = (float) Math.sin(step*0.5);
			float cos_coeff = (float) Math.cos(step*0.5);
			
			rot_quat = new Quaternion(cos_coeff,
					sin_coeff*rot_vector[0], 
					sin_coeff*rot_vector[1],
					sin_coeff*rot_vector[2]);
			
			rot_quat_conj = rot_quat.getConjQuat();
			
			// the "correction"-vector
			corr = new float[]{c1.getX(),c1.getY(),c1.getZ()};
			
			
			if (bond_index!=0) {
				this.just_used = new LinkedList<Atom>(just_tested);
				just_used.addAll(disj.get(bonds[bond_index-1]));
				}
			else {
				this.just_used = just_tested;
			}
			
			while (n<no_steps) {
				
				// more to rotate or last bond ?
				if (bond_index!=bonds.length-1) 
					new RotationNode(this,molecule,bonds,rot_map,disj,bond_index+1,n,just_used,fragment_score);
				else 
					new TerminalNode(this,molecule,n,disj.get(bonds[bond_index]),just_used,fragment_score,bond_index);
				
				// rotate 
				for (Atom a:rot_map.get(bonds[bond_index]))
					a.rotate(rot_quat,rot_quat_conj, corr);
				
				n++;
			}
			
			this.just_used = null;
			
			recomputeBorders(false);
			
		}
		
		public void rotate(int steps,Map<Bond,List<Atom>> rot_map) {
			
			int n = 0;
			while (n<steps) {
				for (Atom a:rot_map.get(bonds[bond_index]))
					a.rotate(rot_quat,rot_quat_conj, corr);
				n++;
			}
		}
	}
	
	
	class TerminalNode extends ConfSearchTree {
		
		protected final int step_before;
		
		public TerminalNode(ConfSearchTree parent,Molecule molecule,int step_before,Set<Atom> last_atoms,List<Atom> just_tested,double[] fragment_score,int bond_index) {
			super(parent,molecule); // new Molecule() leads to enourmous GC
			this.step_before = step_before;
			
			double temp = 0;
			if (last_atoms==null) {
				min_score = molecule.getFastScore(E);
			}
			else {
				temp = getScore(E,molecule,last_atoms,just_tested);
				if (Double.isNaN(temp)) {
					prune();
					return;
				}
				min_score = parent.part_score + temp + fragment_score[bond_index];
			}
			
			max_score = min_score;
		}
		
		public int size() {
			return 1;
		}
		
	}
	
}
