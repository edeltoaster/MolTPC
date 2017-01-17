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
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class TautIt implements Iterable<Molecule> {

	private final Molecule root;
	private final TautRegion[] tautRegions;
	private final TrTries tries;
	private final Map<Atom,Set<Atom>> altReachMap;
	private final Map<DonAccPair,List<Integer>> protonMove;
	private final boolean aromatic_breakage;
	private final boolean sp2_preservation;
	private final boolean stereo;
	private boolean has_aromatic = false;
	List<Molecule> molecules = null;
	List<TautIt> opened_chains = new LinkedList<TautIt>();
	
	
	public TautIt(Molecule root,boolean aromatic_breakage,boolean sp2_preservation,boolean stereo,boolean chain_opening,boolean verbose) {
		this.root = root;
		this.aromatic_breakage = aromatic_breakage;
		this.sp2_preservation = sp2_preservation;
		this.stereo = stereo;
		for (Atom a:root.getFrame()) {
			if (a.isAromatic()) {
				has_aromatic = true;
				break;
			}
		}
	
		this.altReachMap = root.getAltReachMap(aromatic_breakage,true);
		this.tautRegions = preprocess();
		
		// debugging outs
		if (verbose) {
			for (int i=0;i<tautRegions.length;i++)
				System.out.println("Site "+(i+1)+":"+tautRegions[i]);
			printDebugOut();
		}
		
		this.protonMove = preprocessBonds();
		this.tries = buildTries();
		
		if (chain_opening) {
			for (Molecule m:this.detectRingChain())
				opened_chains.add(new TautIt(m, aromatic_breakage, sp2_preservation, stereo, false,false));
			if (verbose && opened_chains.size()>0)
				System.out.println(opened_chains.size()+" variants with open rings.");
		}
	}
	
	/**
	 * searches, identifies and saves all "tautomeric regions"
	 * by merging the altReachability-Mapping from Molecule;
	 * as a side-effect, bonds used for proton-movement are saved
	 */
	private TautRegion[] preprocess() {
		List<TautRegion> regions = new LinkedList<TautRegion>();
		Map<Atom,Set<Atom>> map = root.getFilteredAltReachMap(aromatic_breakage,sp2_preservation);
		Atom[] don = map.keySet().toArray(new Atom[map.keySet().size()]);
		List<Atom> used = new LinkedList<Atom>();
		
		//temp-lists used in loops
		List<Atom> don_temp;
		List<Atom> acc_temp;
		Set<Atom> shared_temp;
		
		//needed for persistence
		Set<Atom> shared_persistent = new HashSet<Atom>();
		
		//try to merge donors in same region
		for (int i=0;i<don.length;i++) {
			
			//donor was merged into another region before
			if (used.contains(don[i]))
				continue;
			
			don_temp = new LinkedList<Atom>();
			don_temp.add(don[i]);
			acc_temp = new LinkedList<Atom>(map.get(don[i]));
			shared_temp = new HashSet<Atom>();
			
			for (int j=i+1;j<don.length;j++) {
				
				//check: totally mergable?
				if ((map.get(don[i]).containsAll(map.get(don[j]))&&
						map.get(don[j]).containsAll(map.get(don[i])))) {
					don_temp.add(don[j]);
					used.add(don[j]);
					//System.out.println("Merge!");
					continue;
				}
				
				//computing intersection
				Set<Atom> temp = new HashSet<Atom>(map.get(don[i]));
				temp.retainAll(map.get(don[j]));
				
				//check: totally disjunct?
				if (temp.size()==0) {
					//System.out.println("Disjunct!");
					continue;
				}
				
				//add shared endpoints
				//System.out.println("shared points!");
				shared_persistent.addAll(temp);
			}
			
			//shared endpoints check
			for (Atom a:acc_temp)
				if (shared_persistent.contains(a))
					shared_temp.add(a);
					
			//add processed region
			regions.add(new TautRegion(don_temp.toArray(new Atom[don_temp.size()]),
					acc_temp.toArray(new Atom[acc_temp.size()]),
					shared_temp.toArray(new Atom[shared_temp.size()])));
			
		}
		
		map = null;
		
		return regions.toArray(new TautRegion[regions.size()]);
	}
	
	private Map<DonAccPair,List<Integer>> preprocessBonds() {
		Map<DonAccPair,List<Integer>> map = new HashMap<TautIt.DonAccPair, List<Integer>>();
		
		DonAccPair pair;
		
		for (TautRegion t:tautRegions) 
			for (Atom d:t.getDonors())
				for (Atom a:t.getAcceptors()) {
					pair = new DonAccPair(d.getNo(),a.getNo());
					map.put(pair,getShBondPath(d,a));
					root.unMarkAtoms();
				}
		
		return map;
	}
	
	private TrTries buildTries() {
		TrTries root = new TrTries(this.root);
		
		if (tautRegions.length==0)
			return root;
		
		// action from hell
		List<TrTries> temp = new LinkedList<TautIt.TrTries>();
		for (TautRegion r:tautRegions) {
			root.getNodes(temp);
			for (TrTries tr:temp)
				r.buildAndAddSubTries(tr);
			temp.clear();
		}
		
		return root;
	}
	
	public List<Molecule> detectRingChain() {
		List<Molecule> opened = new LinkedList<Molecule>();
		List<int[]> variants = new LinkedList<int[]>();
		
		for (Atom a:root.getCycleAtoms()) {
			if (a.getType() == Element.O) {
				for (Atom b:a.getAdjFrame(1,false))
					if (b.getType() == Element.C)
						for (Atom c:b.getAdjFrame(1,false)) {
							if (c == a)
								continue;
							if (c.isOH())
								variants.add(new int[]{a.getNo(),b.getNo(),c.getNo()});
						}
			}
		}
		
		List<Molecule> temp;
		opened.add(root);
		for (int[] i:variants) {
			temp = new LinkedList<Molecule>();
			
			for (Molecule m:opened)
				temp.add(new Molecule(m,i[0],i[1],i[2]));
			
			opened.addAll(temp);
		}
		
		opened.remove(0);
		return BabelInterface.batch_sd_minimize(opened, 100, 1);
	}
	
	public void printDebugOut() {
		System.out.println("alt-reach:");
		boolean first = true;
		
		for (Atom a:this.altReachMap.keySet()) {
			String temp = a.getNo()+":";
			for (Atom b:this.altReachMap.get(a)) {
				if (first) {
					temp += b.getNo();
					first = false;
				}
				else
					temp += ","+b.getNo();
			}
			System.out.println(temp);
		}
		
		System.out.println("alt-reach filtered:");
		Map<Atom,Set<Atom>> filtered = this.root.getFilteredAltReachMap(aromatic_breakage, sp2_preservation);
		for (Atom a:filtered.keySet()) {
			String temp = a.getNo()+":";
			first = true;
			for (Atom b:filtered.get(a)) {
				if (first) {
					temp += b.getNo();
					first = false;
				}
				else
					temp += ","+b.getNo();
			}
			System.out.println(temp);
		}
	}
	
	/**
	 * represents a region of tautomeric forms
	 *
	 */
	class TautRegion {
		private final Atom[] donors, acceptors, shared;
		
		public TautRegion(Atom[] donors,Atom[] acceptors,Atom[] shared) {
			this.donors = donors;
			this.acceptors = acceptors;
			this.shared = shared;
		}
		
		public TrTries buildAndAddSubTries(TrTries root) {
			List<TrTries> children = new LinkedList<TautIt.TrTries>();
			
			//nulltransfer
			if (donors.length>1) {
				Atom[] new_don = Arrays.copyOfRange(donors,1,donors.length);
				Atom[] new_acc = Arrays.copyOfRange(acceptors,0,acceptors.length);
				TautRegion new_taut = new TautRegion(new_don, new_acc, shared);
				new_taut.buildAndAddSubTries(root);
			}
			
			boolean dead_end,add;
			
			for (int i=0;i<acceptors.length;i++) {
				
				dead_end = false;
				add = false;
				
				//shared-check!
				for (Atom sh:shared)
					if (acceptors[i].equals(sh)) {
						if (root.usedShared(sh))
							dead_end = true;
						else
							add = true;
					}
				
				// dead end, shared acceptor already used
				if (dead_end)
					continue;
				
				TrTries temp = new TrTries(root.getMolecule(),root,donors[0].getNo(),acceptors[i].getNo());
				
				// cut if it makes no sense going further for plausible molecules
				               
					if (add)
						temp.addShared(acceptors[i]);
				
					if (donors.length>1 && ((i+1)<acceptors.length)) {
						Atom[] new_don = Arrays.copyOfRange(donors,1,donors.length);
						Atom[] new_acc = Arrays.copyOfRange(acceptors,i+1,acceptors.length);
						TautRegion new_taut = new TautRegion(new_don, new_acc, shared);
						new_taut.buildAndAddSubTries(temp);
					}
				
					children.add(temp);
					
					// stereovariant possible??
					if (!stereo)
						continue;
					
					boolean is_stereo = false;
					for (Atom a:temp.getMolecule().getChiralAtoms()) {
						if (a.getNo() == acceptors[i].getNo())
							is_stereo = true;
					}
					
					// if yes:
					if (is_stereo) {
						temp = new TrTries(root.getMolecule(),root,donors[0].getNo(),acceptors[i].getNo());
						Molecule m = temp.getMolecule();
						m.inverseH(m.getAtomNo(acceptors[i].getNo()).getAdjHs().get(0),m.getAtomNo(acceptors[i].getNo()));
						                  
						if (add)
							temp.addShared(acceptors[i]);
					
						if (donors.length>1 && ((i+1)<acceptors.length)) {
							Atom[] new_don = Arrays.copyOfRange(donors,1,donors.length);
							Atom[] new_acc = Arrays.copyOfRange(acceptors,i+1,acceptors.length);
							TautRegion new_taut = new TautRegion(new_don, new_acc, shared);
							new_taut.buildAndAddSubTries(temp);
						}
					
						children.add(temp);
					}
				
				
			}
			
			root.getChildren().addAll(children);
			return root;
		}
		
		public Atom[] getDonors() {
			return donors;
		}
		
		public Atom[] getAcceptors() {
			return acceptors;
		}
		
		public Atom[] getShared() {
			return shared;
		}
		
		public String toString() {
			String temp = Integer.toString(donors[0].getNo());
			
			for (int i=1;i<donors.length;i++)
				temp += ","+donors[i].getNo();
			
			temp += " over ";
			
			temp += Integer.toString(acceptors[0].getNo());
			
			for (int i=1;i<acceptors.length;i++)
				temp += ","+acceptors[i].getNo();
			
			if (shared.length>0) {
				temp += " shared:"+shared[0].getNo();
				for (int i=1;i<shared.length;i++)
					temp += ","+shared[i].getNo();
			}
			return temp;
		}
		
	}
	
	class TrTries {
		private TrTries parent;
		final private List<TrTries> children = new LinkedList<TautIt.TrTries>();
		final private Set<Atom> used_shared;
		final private Molecule molecule;
		
		// root constructor
		public TrTries(Molecule mol) {
			this.molecule = mol;
			this.used_shared = new HashSet<Atom>();
			this.parent = null;
		}
		
		// non-root constructor
		public TrTries(Molecule mol,TrTries parent,int don,int acc) {
			this.molecule = new Molecule(mol);
			this.parent = parent;
			this.molecule.transform(don, acc,protonMove.get(new DonAccPair(don,acc)));
			this.used_shared = new HashSet<Atom>(this.parent.used_shared);
			this.molecule.setName(mol.getName()+"("+don+"->"+acc+")");
		}
		
		public List<TrTries> getNodes(List<TrTries> temp) {
			temp.add(this);
			for (TrTries t:children)
				t.getNodes(temp);
			
			return temp;
		}
		
		public List<Molecule> getMolecules(List<Molecule> molecules) {
			molecules.add(this.molecule);
			for (TrTries t:children)
				t.getMolecules(molecules);
			
			return molecules;
		}
		
		public int size() {
			int sum = 1;
			for (TrTries tr:children)
				sum += tr.size();
			
			return sum;
		}
		
		public List<TrTries> getChildren() {
			return children;
		}
		
		public TrTries getParent() {
			return parent;
		}
		
		public void setParent(TrTries parent) {
			this.parent = parent;
		}
		
		public Molecule getMolecule() {
			return molecule;
		}
		
		public void addShared(Atom sh) {
			used_shared.add(sh);
		}
		
		public boolean usedShared(Atom a) {
			return used_shared.contains(a);
		}
		
	}
	
	class DonAccPair {
		final private int donor,acceptor;
		
		public DonAccPair(int donor,int acceptor) {
			this.donor = donor;
			this.acceptor = acceptor;
		}

		public int getDonor() {
			return donor;
		}

		public int getAcceptor() {
			return acceptor;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + acceptor;
			result = prime * result + donor;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (!(obj instanceof DonAccPair))
				return false;
			DonAccPair other = (DonAccPair) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (acceptor != other.acceptor)
				return false;
			if (donor != other.donor)
				return false;
			return true;
		}

		private TautIt getOuterType() {
			return TautIt.this;
		}
		
	}
	
	public Molecule getRoot() {
		return root;
	}

	TautRegion[] getTautRegions() {
		return tautRegions;
	}
	/**
	 * calculates the number of tautomeric forms
	 * @return number of forms
	 */
	public int size() {
		return tries.size();
	}
	
	public List<Molecule> getAllForms() {
		if (molecules != null)
			return molecules;
		
		List<Molecule> temp_molecules = new LinkedList<Molecule>();
		molecules = new LinkedList<Molecule>();
		tries.getMolecules(temp_molecules);
		
		// get opened-chain variants
		for (TautIt ti:opened_chains)
			ti.tries.getMolecules(temp_molecules);
		
		for (Molecule m:temp_molecules) {
			if (aromatic_breakage && has_aromatic) {
				m.resetAromatic();
				if (m.QaDValencyFix())
					molecules.add(m);
			} else
				molecules.add(m);	
		}
		
		Set<String> seen_molecules = new HashSet<String>();
		List<Molecule> pre_output = BabelInterface.batch_sd_minimize(molecules, 100, 1);
		
		molecules.clear();
		
		// batch much faster ...
		List<String> smi = BabelInterface.batch_SMILES(pre_output);
		Iterator<String> smi_it = smi.iterator();
		Iterator<Molecule> mol_it = pre_output.iterator();
		
		Molecule m;
		String m_smi;
		while (mol_it.hasNext()) {
			m = mol_it.next();
			m_smi = smi_it.next();
			
			if (seen_molecules.contains(m_smi))
				continue;
			
			molecules.add(m);
			seen_molecules.add(m_smi);
			
		}
		
		return molecules;
	}
	
	private List<Integer> getShBondPath(Atom d,Atom a) {
		Deque<Atom> path = new LinkedList<Atom>();
		Atom temp = d;
		List<Atom> next_list = new LinkedList<Atom>();
		
		//loop
		while (!temp.equals(a)) {
			path.push(temp);
			temp.mark();
			next_list.clear();
			
			for (Atom next:temp.getAdjFrame(0,true)) {
				if (next.isMarked() || !altReachMap.get(d).contains(next))
					continue;
				next_list.add(next);
			}
			
			if (next_list.size()==0) {
				path.pop();
				temp = path.pop();
			}
			else {
				temp = next_list.get(0);
			}
		}
		// loop end
		
		path.push(a);
		
		List<Integer> list = new LinkedList<Integer>();
		for (Atom p:path)
			list.add(p.getNo());
		return list;
	}
	
	@Override
	public Iterator<Molecule> iterator() {
		
		return getAllForms().iterator();
	}

}
