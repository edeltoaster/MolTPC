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
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;


public class EnsembleAnalyzer {
	public static final double R = 8.314472;
	private Map<Double,Molecule> ensemble;

	public EnsembleAnalyzer(Collection<Molecule> ensemble,float epsilon,int first_min,int no_best,int sec_min) {
		min_and_trim(ensemble,epsilon,first_min,no_best, sec_min);
	}
	
	public EnsembleAnalyzer(Collection<Molecule> ensemble,float epsilon) {
		this.ensemble = BabelInterface.batch_energy_sorted(ensemble, epsilon);
	}
	
	public EnsembleAnalyzer(Collection<Molecule> ensemble,float epsilon,boolean cg) {
		if (cg) { // conjugate gradient
			this.ensemble = BabelInterface.batch_energy_sorted(BabelInterface.batch_cg_minimize(ensemble,15000,epsilon), epsilon);
		} else { // steepest descent
			this.ensemble = BabelInterface.batch_energy_sorted(BabelInterface.batch_sd_minimize(ensemble,15000,epsilon), epsilon);
		}
	}

	private void min_and_trim(Collection<Molecule> ensemble,float epsilon,int first_min,int no_best,int sec_min) {
		Map<Double,Molecule> first = BabelInterface.batch_energy_sorted(BabelInterface.batch_sd_minimize(ensemble,first_min,epsilon), epsilon);
		
		// getting best ones for second minimizing step
		List<Molecule> temp = new LinkedList<Molecule>();
		int n = 0;
		for (double f:first.keySet()) {
			if (n>=no_best)
				break;
			temp.add(first.get(f));
			n++;
		}
		
		first = null;
		this.ensemble = BabelInterface.batch_energy_sorted(BabelInterface.batch_cg_minimize(temp,sec_min,epsilon), epsilon);
	}
	
	public void cutBoltzmann(int partition,int temperature) {
		if (partition<=0)
			throw new IllegalArgumentException("cutBoltzmann wrong pratition-arg");
		
		float max_dE = (float) (- Math.log(1/(double)partition) * temperature * R /4184.0); // in kcal / mol
		cut(max_dE);
	}
	
	public void cutThermic(int temperature) {
		float max_dE = (float) (3 * temperature * R /4184.0); // in kcal/mol
		cut(max_dE);
	}
	
	public void cut(float max_dE) {
		double min_energy = Double.MAX_VALUE;
		Map<Double,Molecule> temp = new TreeMap<Double, Molecule>();
		
		for (double f:ensemble.keySet()) {
			if (f<min_energy) // only happens for the first float, since datastructure is a SORTED treeset
				min_energy = f;
			
			// since here -> too bad
			if (f>(min_energy+max_dE))
				break;
			
			temp.put(f,ensemble.get(f));
		}
		
		ensemble = temp;
	}
	
	public Map<Double,Molecule> getEnsemble() {
		return ensemble;
	}
	
	public int size() {
		return ensemble.size();
	}
	
	public void writeEnsembleToFolder(String folder) {
		int n = 0;
		String temp = null;
		
		if (folder.endsWith("/"))
			temp = folder;
		else
			temp = folder+"/";
		
		for (double f:ensemble.keySet()) {
			ensemble.get(f).writeHinFile(temp+(n++)+"-"+f+".hin");
		}
	}
	
	public void printConfigurations() {
		Map<String,Double> map = getConfigurations();
		System.out.println("Distribution of configurations:");
		for (String s:map.keySet()) {
			System.out.printf("%5.2f%% : %s%n",map.get(s),s);
		}
	}
	
	public Map<String,Double> getConfigurations() {
		String smiles;
		Map<String,Double> map = new HashMap<String,Double>();
		for (Molecule m:ensemble.values()) {
			smiles = BabelInterface.getSMILESfast(m);
			if (map.containsKey(smiles))
				map.put(smiles,map.get(smiles)+1);
			else
				map.put(smiles,1.0);
		}
		
		for (String s:map.keySet()) {
			map.put(s,map.get(s)/ensemble.size()*100);
		}
		
		BabelInterface.smiClean();
		
		return map;
	}
	
	public void printConfigurationsBoltzmann(double T) {
		
		Map<String,Double> map = getConfigurationsBoltzmann(T);
		
		System.out.println("Distribution of configurations (Boltzmann-weighted):");
		for (String s:map.keySet()) {
			System.out.printf("%5.2f%% : %s%n",map.get(s),s);
		}
		
	}
	
	public Map<String,Double> getConfigurationsBoltzmann(double T) {
		String smiles;
		double sum = 0;
		double temp = 0;
		double min = Double.POSITIVE_INFINITY;
		
		//factor
		double factor = R/4184.0 * T; // "boltzmann-factor"
		Map<String,Double> map = new HashMap<String,Double>();
		
		for (double e:ensemble.keySet()) {
			// set min
			if (Double.isInfinite(min)) {
				min = e;
			}
			
			smiles = BabelInterface.getSMILESfast(ensemble.get(e));
			
			// bla boltzmann:
			temp = 1/Math.exp((e-min)/factor);
			sum += temp;
			if (map.containsKey(smiles))
				map.put(smiles,map.get(smiles)+temp);
			else
				map.put(smiles,temp);
			
		}
		
		for (String s:map.keySet()) {
			map.put(s,map.get(s)/sum*100);
		}
		
		BabelInterface.smiClean();
		
		return map;
	}
	
	public static void printConfClusters(double T,int threshold,List<Map<Double,Molecule>> ensembles,String repr_path) {
		
		double factor = R/4184.0 * T; // "boltzmann-factor"
		int no_runs = ensembles.size();
		
		System.out.println("Distribution of conformations:");
		
		// sort by configurations
		Map<String,List<Molecule>> inchi_mol_map = new HashMap<String,List<Molecule>>();
		Map<Molecule,Double> mol_contr_map = new HashMap<Molecule, Double>();
		String inchi;
		Molecule temp_mol;
		double sum,temp_contr,min;
		
		for (Map<Double,Molecule> ens:ensembles) {
			min = Double.POSITIVE_INFINITY;
			sum = 0;
			for (double e:ens.keySet()) {
				// set min
				if (Double.isInfinite(min)) {
					min = e;
				}
				temp_mol = ens.get(e);
				inchi = BabelInterface.getINCHIfast(temp_mol);
				
				// compute contribution
				temp_contr = 1/Math.exp((e-min)/factor);
				sum += temp_contr;
				mol_contr_map.put(temp_mol,temp_contr);
				if (inchi_mol_map.containsKey(inchi))
					inchi_mol_map.get(inchi).add(temp_mol);
				else {
					List<Molecule> temp_list = new LinkedList<Molecule>();
					temp_list.add(temp_mol);
					inchi_mol_map.put(inchi,temp_list);
				}
			}
			
			// normalizing
			for (Molecule m:ens.values()) {
				mol_contr_map.put(m,(mol_contr_map.get(m)/sum/no_runs*100));
			}
		}
		// molecules now separated by INCHI and contribution per ensemble in mol_contr_map
		
		// cluster conformations of configurations
		double out_sum = 0;
				for (String inch:inchi_mol_map.keySet()) {
					
					// inits
					List<List<Float[]>> cluster = new LinkedList<List<Float[]>>();
					System.out.println(BabelInterface.getSMILESfast(inchi_mol_map.get(inch).get(0)));
					List<Float[]> samples = new LinkedList<Float[]>();
					Map<Float[],Double> angles_contr_map = new HashMap<Float[], Double>();
					Map<Float[],Molecule> angles_to_mol_map = new HashMap<Float[], Molecule>();
					Set<Bond> bonds = null;
					
					// compute angles
					for (Molecule m:inchi_mol_map.get(inch)) {
						bonds = m.getRotatableBondsMap().keySet();
						Float[] angles = new Float[bonds.size()];
						int i = 0;
						for (Bond b:bonds) {
							angles[i++] = (float) b.getAngle();
						}
						samples.add(angles);
						angles_contr_map.put(angles,mol_contr_map.get(m));
						angles_to_mol_map.put(angles,m);
					}
				
					// action
					Deque<Float[]> stack = new LinkedList<Float[]>();
					Float[] current;
					List<Float[]> current_cluster, samples_new;
					while (!samples.isEmpty()) {
						stack.push(samples.remove(0));
						current_cluster = new LinkedList<Float[]>();
						
						while (!stack.isEmpty()) {
							current = stack.pop();
							current_cluster.add(current);
							samples_new = new LinkedList<Float[]>();
							
							for (Float[] f:samples) {
								if (areSimiliar(current, f, threshold)) 
									stack.push(f);
								else
									samples_new.add(f);
							}
							
							samples = samples_new;
						}
						
						cluster.add(current_cluster);
					}
					
					// output information per cluster per molecule
					System.out.println(cluster.size()+" distinguishable conformation(s) found in "+inchi_mol_map.get(inch).size()+" molecules.");
					int dim = cluster.get(0).get(0).length;
					String temp = "";
					for (Bond b:bonds) {
						int[] atoms = b.getAngleAtoms();
						temp += "["+atoms[0]+"-"+atoms[1]+"-"+atoms[2]+"-"+atoms[3]+" "+"]";
					}
					System.out.println(temp);
					
					Map<List<Float[]>,Double> cl_contr_map = new HashMap<List<Float[]>, Double>();
					List<List<Float[]>> sorted_list = new LinkedList<List<Float[]>>();
					
					for (List<Float[]> l:cluster) {
						// total contribution
						sum = 0;
						for (Float[] f:l) {
							sum += angles_contr_map.get(f);
						}
						cl_contr_map.put(l,sum);
					}
					
				
					// sort
					double max = Double.NEGATIVE_INFINITY;
					List<Float[]> max_element = null;
					while (!cluster.isEmpty()) {
						for (List<Float[]> l:cluster) {
							if (cl_contr_map.get(l)>max) {
								max = cl_contr_map.get(l);
								max_element = l;
							}
						}
						max = Double.NEGATIVE_INFINITY;
						sorted_list.add(max_element);
						cluster.remove(max_element);
					}
					
					// output analysis (until some threshold is reached)
					boolean minor = sorted_list.size() < 8;
					for (List<Float[]> l:sorted_list) {
						// check threshold criteria
						if (!minor && (out_sum > 80 || cl_contr_map.get(l) <= 2))
							break;
						
						out_sum += cl_contr_map.get(l);
						
						// search highest contributing representant
						Float[] highest_contr_repr = null;
						double highest_contr = Double.NEGATIVE_INFINITY;
						for (Float[] angles:l) {
							double c_temp = angles_contr_map.get(angles);
							if (c_temp > highest_contr) {
								highest_contr = c_temp;
								highest_contr_repr = angles;
							}
						}
						
						// avg
						Float[] angles = averagedAngles(l);
						temp = "";
						for (int i=0;i<dim;i++) {
							if (i == dim-1)
								temp += String.format("%6.2f",angles[i]);
							else
								temp += String.format("%6.2f ",angles[i]);
						}
						
						System.out.printf("%5.2f%% : %s (%dx)",cl_contr_map.get(l),temp,l.size());
						
						// makes no sense if there is only one member
						boolean print_more = l.size() > 1;
						
						if (print_more) {
							// compute max deviation per angle
							angles = maxAngleDeviationRANGE(l);
							print_more = false; 
							for (Float dev:angles) // check if molecules have deviation
								if (dev != 0.0) {
									print_more = true;
									break;
								}	
						}
						
						// printing deviations makes sense
						if (print_more) {
							// from range
							temp = "";
							for (int i=0;i<dim;i++) {
								if (i == dim-1)
									temp += String.format("%5.2f",angles[i]);
								else
									temp += String.format("%5.2f ",angles[i]);
							}
							System.out.printf(" range:%s",temp);
							
							//from avg
							angles = maxAngleDeviationAVG(l);
							temp = "";
							for (int i=0;i<dim;i++) {
								if (i == dim-1)
									temp += String.format("%6.2f",angles[i]);
								else
									temp += String.format("%6.2f ",angles[i]);
							}
							System.out.printf(" dAVG:%s%n",temp);
							
						} else {
							System.out.printf("%n");
						}
						
						// writeout representative molecule with highest contribution
						String outputfile = String.format("%s%5.2f.hin",repr_path,cl_contr_map.get(l));
						angles_to_mol_map.get(highest_contr_repr).writeHinFile(outputfile);
						
					}
					System.out.println();
				}
				
				BabelInterface.smiClean();
	}
	
	private static boolean areSimiliar(Float[] angles1, Float[] angles2, int threshold) {
		if (threshold <= 0)
			throw new IllegalArgumentException("Threshold for ConfClustering must have a positive value");
		
		if(angles1.length != angles2.length)
			return false;
		
		for (int i=0;i<angles1.length;i++) {
			if (Math.abs(angles1[i]-angles2[i]) > threshold)
				return false;
		}
		
		return true;
	}
	
	private static Float[] averagedAngles(List<Float[]> angles) {
		int dim = angles.get(0).length;
		Float[] avg = new Float[dim]; // init of java primitive types: 0.0, but not container
		for (int i=0;i<dim;i++)
			avg[i] = 0.0f;
		
		// summing up
		for (Float[] f:angles) {
			for (int i=0;i<dim;i++)
				avg[i] += f[i];
		}
		
		// ... and divide ...
		for (int i=0;i<dim;i++)
			avg[i] /= angles.size(); // divide by sample-size
		
		return avg;
	}
	
	private static Float[] maxAngleDeviationAVG(List<Float[]> angles) {
		int dim = angles.get(0).length;
		Float[] max_dev = new Float[dim]; // init of java primitive types: 0.0, but not container
		for (int i=0;i<dim;i++)
			max_dev[i] = 0.0f;
		
		// avg
		Float[] avg = averagedAngles(angles);
		
		// summing up
		for (Float[] f:angles) {
			for (int i=0;i<dim;i++) {
				if (Math.abs(f[i]-avg[i]) > Math.abs(max_dev[i]))
					max_dev[i] = f[i]-avg[i];
			}
		}
		
		return max_dev;
	}
	
	private static Float[] maxAngleDeviationRANGE(List<Float[]> angles) {
		int dim = angles.get(0).length;
		Float[] max_dev = new Float[dim]; // init of java primitive types: 0.0, but not container
		Float[] min = new Float[dim];
		Float[] max = new Float[dim];
	
		// checking
		for (int i=0;i<dim;i++) {
			min[i] = Float.MAX_VALUE;
			max[i] = Float.MIN_VALUE;
			for (Float[] f:angles) {
				if (f[i]<min[i])
					min[i] = f[i];
				if (f[i]>max[i])
					max[i] = f[i];
			}
			
			max_dev[i] = max[i] - min[i];
		}
		return max_dev;
	}
	
	public void printConformations() {
		String smiles;
		Map<String,List<Molecule>> map = new HashMap<String,List<Molecule>>();
		for (Molecule m:ensemble.values()) {
			smiles = BabelInterface.getSMILESfast(m);
			if (map.containsKey(smiles))
				map.get(smiles).add(m);
			else {
				List<Molecule> temp_list = new LinkedList<Molecule>();
				temp_list.add(m);
				map.put(smiles,temp_list);
			}
		}
		
		System.out.println("Distribution of conformations:");
		
		for (String smi:map.keySet()) {
			System.out.println("for "+map.get(smi).size()+"x "+smi);
			
			Map<Set<Atom>,ConfTree> frag_node_map = new HashMap<Set<Atom>,ConfTree>(); // stores which node/fragment added this job
			Map<Set<Atom>,Atom[]> frag_torsion_map = new HashMap<Set<Atom>, Atom[]>();
			Deque<Set<Atom>> stack = new LinkedList<Set<Atom>>();
			ConfTree root = null;
			
			for (Molecule m:map.get(smi)) {

				Map<Bond,List<Atom>> rot_map = m.getRotatableBondsMap();
				Atom start = m.getFixAtoms(Molecule.getDisjunctRotAtomsMap(rot_map)).get(0);
				List<Set<Atom>> fragments = m.getFragments();
				Set<Bond> rot_bonds = rot_map.keySet();
				Map<Atom,List<Bond>> atom_to_bond = new HashMap<Atom, List<Bond>>();
				Map<Bond,Boolean> bondMarks = new HashMap<Bond, Boolean>();

				for (Bond b:rot_bonds) {
					if (!atom_to_bond.containsKey(b.getA1()))
						atom_to_bond.put(b.getA1(),new LinkedList<Bond>());
					atom_to_bond.get(b.getA1()).add(b);

					if (!atom_to_bond.containsKey(b.getA2()))
						atom_to_bond.put(b.getA2(),new LinkedList<Bond>());
					atom_to_bond.get(b.getA2()).add(b);
					
					bondMarks.put(b,false);
				}

				for (Set<Atom> fr:fragments) // finding start-fragment (the one considered as fixed by confsearch)
					if (fr.contains(start)) {
						stack.add(fr);
						frag_node_map.put(fr,null);
						break;
					}

				// building tree
				Set<Atom> current = null;
				ConfTree current_node = null;
				String current_smiles = null;

				while (stack.size()>0) {
					current = stack.pop();
					current_smiles = BabelInterface.getSMILESfast(new Molecule(current));
					
					// only for very first fragment of all
					if (root == null) {
						root = new ConfTree(current_smiles,current);
						current_node = root;
					}
					// first new fixed fragment
					else if (frag_node_map.get(current) == null) {
						current_node = root;
					}
				    else // all other fragments
				    	for (ConfTree tr:frag_node_map.get(current).getChildren()) {
				    		if (tr.getFragment_smiles().equals(current_smiles)) {
				    			if (!frag_torsion_map.get(current)[0].equals(frag_torsion_map.get(tr.getSet())[0]))
				    				continue;
				    			current_node = tr;
				    			tr.visit();
				    			break;
				    		}
				    	}

					// if this is the first fragment like this
					if (current_node == null)
						current_node = new ConfTree(frag_node_map.get(current),current_smiles,current);

					// if this fragment is not the fixed one -> check angles!
					if (current_node.getParent()!=null) {
						Atom[] torsion_atoms = frag_torsion_map.get(current);
						float angle = (float)(Molecule.getTorsionAngle(torsion_atoms[0],torsion_atoms[1],torsion_atoms[2],torsion_atoms[3]));
						current_node.addAngle(angle);
					}

					// look for adjacent fragments
					for (Atom a:current) {
						if (!atom_to_bond.keySet().contains(a))
							continue;

						List<Atom> to_remove = new LinkedList<Atom>();
						to_remove.add(a);
						for (Bond b:atom_to_bond.get(a)) {
							if (bondMarks.get(b))
								continue;
							bondMarks.put(b,true);
							Atom other = b.getOtherAtom(a);

							// find the "other atom" in a fragment
							for (Set<Atom> fr:fragments)
								if (fr.contains(other)) {
									stack.add(fr);
									frag_node_map.put(fr,current_node);
									Atom[] torsion_atoms = new Atom[4];

									// setting atoms from the bond
									torsion_atoms[1] = a;
									torsion_atoms[2] = other;

									// choosing an adjacent from the corresponding sets
									// from current fragment
									torsion_atoms[0] = a.getAdjAtoms().get(0);
									if (torsion_atoms[0].equals(other))
										torsion_atoms[0] = a.getAdjAtoms().get(1);

									// from other fragments
									torsion_atoms[3] = other.getAdjAtoms().get(0);
									if (torsion_atoms[3].equals(a))
										torsion_atoms[3] = other.getAdjAtoms().get(1);

									frag_torsion_map.put(fr,torsion_atoms);
									break;
								}
						}

						// remove from set so this doesnt happen twice
						for (Atom r:to_remove) 
							atom_to_bond.remove(r);


					}

					current_node = null; // for convenience
				}

			}
			// output in some way....
			System.out.println(root);
		}
	}
	
	public void printFirstAppearance() {
		String smiles;
		Set<String> seen = new HashSet<String>();
		int n=0;
		System.out.println("First appearances of configurations:");
		for (Molecule m:ensemble.values()) {
			smiles = BabelInterface.getSMILESfast(m);
			if (!seen.contains(smiles)) {
				seen.add(smiles);
				//output
				System.out.println(n+": "+smiles);
			}
			n++;
		}
		
		BabelInterface.smiClean();
	}
	
	
	class ConfTree {
		private ConfTree parent;
		private List<ConfTree> children = new LinkedList<EnsembleAnalyzer.ConfTree>();
		private List<Float> angles = new LinkedList<Float>();
		private int visits = 1;
		private String fragment_smiles;
		private boolean completed = false;
		private Set<Atom> set;
		
		ConfTree() {
			parent = null;
			fragment_smiles = "";
			this.set = null;
		}
		
		ConfTree(ConfTree parent,String smiles,Set<Atom> set) {
			this.parent = parent;
			if (parent!=null)
				parent.children.add(this);
			this.fragment_smiles = smiles;
			this.set = set;
		}
		
		ConfTree(String smiles,Set<Atom> set) {
			this.parent = null;
			if (parent!=null)
				parent.children.add(this);
			this.fragment_smiles = smiles;
			this.set = set;
		}
		
		public ConfTree getParent() {
			return parent;
		}

		public void setParent(ConfTree parent) {
			this.parent = parent;
		}

		public String getFragment_smiles() {
			return fragment_smiles;
		}
		
		public void setFragment_smiles(String fragment_smiles) {
			this.fragment_smiles = fragment_smiles;
		}

		public List<ConfTree> getChildren() {
			return children;
		}

		public List<Float> getAngles() {
			return angles;
		}
		
		public void addAngle(float angle) {
			angles.add(angle);
		}
		
		public Set<Atom> getSet() {
			return set;
		}
		
		public void visit() {
			visits++;
		}
		
		public int getVisits() {
			return visits;
		}
		
		public void setCompleted() {
			completed = true;
		}
		
		public boolean isCompleted() {
			return completed;
		}
		
		public String toString() {
			String temp = "";
			//for (float angle:angles)
			//	temp += (angle+",");
			temp += fragment_smiles;
			if (children.size()>0){
				temp += "[";
				
				for (ConfTree tr:children)
					temp += tr.toString()+",";
						
				temp += "]";
			}
				
			return temp;
		}
	}
}
