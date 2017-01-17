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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;




public class ProtIt implements Iterable<Molecule> {
	
	private final Molecule root;
	private final int nettocharge;
	private final boolean removeHs;
	private final List<Molecule> protStates;
	private final int[] score;
	private boolean more_to_evaluate = true;
	private boolean[] alternatives;
	private int mol_size;
	
	
	public ProtIt(Molecule molecule,int nettocharge) {
		if (!molecule.isExplHs())
			throw new UnsupportedOperationException("protonation states without explicit protons not possible");
		this.root = molecule;
		molecule.neutralizeCharges();
		this.nettocharge = nettocharge;
		this.removeHs = (nettocharge<0);
		this.score = new int[molecule.getFrame().length];
		this.alternatives = new boolean[score.length];
		this.mol_size = score.length;
		
		if (nettocharge == 0) {
			this.protStates = new LinkedList<Molecule>();
			this.protStates.add(molecule);
		}
		else if (!removeHs)
			this.protStates = BabelInterface.batch_sd_minimize(buildProtonationStates_positive(),100,1);
		else
			this.protStates = BabelInterface.batch_sd_minimize(buildProtonationStates_negative(),100,1);
		
	}
	
	private void initRun_negative(boolean reset_alt) {
		Atom[] atoms = root.getFrame();
	
		// init scores[], reset alternatives if needed
		for (int i=0;i<score.length;i++) {
			if (reset_alt)
				alternatives[i] = false;
			if (atoms[i].getAdjHs().size()==0)
				score[i] = Integer.MAX_VALUE;
			else if (atoms[i].isCOOH())
				score[i] = -2 * Math.abs(nettocharge) * mol_size;
			else
				score[i] = atoms[i].getType().getLoseHEnergy() * Math.abs(nettocharge) * mol_size;
		}
	}
	
	private void initRun_positive(boolean reset_alt) {
		Atom[] atoms = root.getFrame();
	
		// init scores[], reset alternatives if needed
		for (int i=0;i<score.length;i++) {
			if (reset_alt)
				alternatives[i] = false;
		if (!atoms[i].onlySB())
			score[i] = Integer.MIN_VALUE;
		else if (atoms[i].isNH2X())
			score[i] = 2 * Math.abs(nettocharge) * mol_size;
		else
			score[i] = atoms[i].getType().getLoseHEnergy() * Math.abs(nettocharge) * mol_size;
		}
	}
	
	private List<Integer> singleRun(int start_index) {
		more_to_evaluate = false;
		List<Integer> Hs = new ArrayList<Integer>(Math.abs(nettocharge)); // for faster equality-checking later, much better in this case than Molecule.equals
		int n = 0;
		
		while (n != Math.abs(nettocharge)) {
			
			Hs.add(start_index);
			
			n++;
			if(!removeHs)
				score[start_index] = Integer.MIN_VALUE;
			else
				score[start_index] = Integer.MAX_VALUE;
			
			// early break
			if (n == Math.abs(nettocharge))
				break;
			
			// update scores of neighbourhood
			Deque<Atom> stack1 = new LinkedList<Atom>();
			Deque<Atom> stack2 = new LinkedList<Atom>();
			stack1.push(root.getFrame()[start_index]);
			
			int w = 1;
			Atom temp;
			
			while (stack1.size()>0 && w<mol_size) {
				
				// new adjacent atoms w bonds away from atom just chosen
				while (stack1.size()>0) {
					temp = stack1.pop();
					for (Atom b:temp.getAdjFrame(0,true))
						if (!b.isMarked()) {
							b.mark();
							stack2.add(b);
						}
				}
				
				// modify score[]
				while (stack2.size()>0) {
					int d_score = (mol_size-w);
					temp = stack2.pop();
					
					if (!removeHs) {
						if (score[temp.getNo()-1] != Integer.MIN_VALUE)
							score[temp.getNo()-1] -= d_score;
					} else {
						if (score[temp.getNo()-1] != Integer.MAX_VALUE)
							score[temp.getNo()-1] += d_score;
					}
					
					stack1.push(temp);
				}
				
				w++;
			}
			
			root.unMarkAtoms();
			
			// choose next atom 
			start_index = 0;
			if (!removeHs) {
				for (int i=1;i<score.length;i++) {
					 // choose lowest score
					if (score[i]>score[start_index]) {
						start_index = i;
					}
					
					// handling alternatives
					if (score[i]==score[start_index]) {
						if (alternatives[i]) {
							start_index = i;
						}
						else {
							alternatives[i] = true;
							more_to_evaluate = true;
						}
					}
				}
			} else {
				for (int i=1;i<score.length;i++) {
					 // choose lowest score
					if (score[i]<score[start_index]) {
						start_index = i;
					}
					
					// handling alternatives
					if (score[i]==score[start_index]) {
						if (alternatives[i]) {
							start_index = i;
						}
						else {
							alternatives[i] = true;
							more_to_evaluate = true;
						}
					}
				}
			}
			
		}
		
		Collections.sort(Hs);
		return Hs;
	}
	
	private List<Molecule> buildProtonationStates_negative() {
		
		Set<List<Integer>> set_of_Hs_to_del = new HashSet<List<Integer>>();
		
		int[] no_first = new int[score.length];
		
		while (true) {
			
			initRun_negative(true); // resets scores and alternatives
			root.unMarkAtoms();
			
			// choosing start_index depending on score and no_first
			int start_index = 0;
			for (int i=1;i<score.length;i++) {
				if (score[i]<score[start_index]) {
				 start_index = i;
				}
				else if (score[i]==score[start_index]) {
					if (no_first[i]<no_first[start_index]) {
						start_index = i;
					}
				}
			}
		
			no_first[start_index]++;
		
			// if no_first[start_index] == 2 -> terminate, probably dependancy on nettocharge?
			if (no_first[start_index]==2)
				break;
			 
			more_to_evaluate = true; // for first run
			while (more_to_evaluate) {// more_to_evaluate is set to false in singleRun and reset to true if there are branches
				set_of_Hs_to_del.add(singleRun(start_index));
				initRun_negative(false); // only reset scores
				}
			 
		 }
		
		// building distinct molecules
		List<Molecule> molecules = new LinkedList<Molecule>();
		
		for (List<Integer> b:set_of_Hs_to_del) {
			Molecule mol = new Molecule(root);
			for (int i:b) {
				Atom losingH = mol.getFrame()[i];
				Atom H = losingH.getAdjHs().get(0);
				losingH.setCharge(-1);
				losingH.resetAtomType();
				mol.deleteH(H);
				
			}
			molecules.add(mol);
		}
		
		
		return molecules;
	}
	
	private List<Molecule> buildProtonationStates_positive() {
		
		Set<List<Integer>> set_of_Hs_to_add = new HashSet<List<Integer>>();
		
		int[] no_first = new int[score.length];
		
		while (true) {
			
			initRun_positive(true); // resets scores and alternatives
			root.unMarkAtoms();
			
			// choosing start_index depending on score and no_first
			int start_index = 0;
			for (int i=1;i<score.length;i++) {
				if (score[i]>score[start_index]) {
				 start_index = i;
				}
				else if (score[i]==score[start_index]) {
					if (no_first[i]<no_first[start_index]) {
						start_index = i;
					}
				}
			}
		
			no_first[start_index]++;
		
			// if no_first[start_index] == 2 -> terminate, probably dependancy on nettocharge?
			if (no_first[start_index]==2)
				break;
			 
			more_to_evaluate = true; // for first run
			while (more_to_evaluate) {// more_to_evaluate is set to false in singleRun and reset to true if there are branches
				set_of_Hs_to_add.add(singleRun(start_index));
				initRun_positive(false); // only reset scores
				}
			 
		 }
		
		// building distinct molecules
		List<Molecule> molecules = new LinkedList<Molecule>();
		
		for (List<Integer> b:set_of_Hs_to_add) {
			Molecule mol = new Molecule(root);
			for (int i:b) {
				Atom gettingH = mol.getFrame()[i];
				mol.addNewH(gettingH);
				gettingH.setCharge(1);
				gettingH.resetAtomType();
				
			}
			molecules.add(mol);
		}
		
		
		return molecules;
	}

	public int size() {
		return this.protStates.size();
	}
	
	public List<Molecule> getAllForms() {
		return protStates;
	}
	
	@Override
	public Iterator<Molecule> iterator() {
		return this.protStates.iterator();
	}

}
