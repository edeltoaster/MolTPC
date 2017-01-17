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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;


public class Molecule implements Comparable<Molecule>{
	private Atom[] frame;
	private List<Atom> protons;
	private String name;
	private boolean explHs;
	private Set<Bond> flipped;
	private Map<Atom,Atom> from;
	private Set<Atom> cycle_atoms;
	public Map<Atom,List<Atom>> score_map;
	public Map<Atom,List<Atom>> score_map_4;
	private List<Atom> chiral_atoms;
	
	public Molecule(Molecule mol) {
		//getting ready
		this.frame = new Atom[mol.frame.length];
		this.protons = new LinkedList<Atom>();
		this.name = new String(mol.name);
		this.explHs = mol.explHs;
		
		//filling with NEW (!) objects
		for (int i=0;i<mol.frame.length;i++) 
			this.frame[i] = new Atom(mol.frame[i]);
		for (Atom a:mol.getProtons())
			this.protons.add(new Atom(a));
		
		//adding NEW (!) bonds
		Set<Bond> bonds = new HashSet<Bond>();
		for (int i=0;i<mol.frame.length;i++) 
			for (Bond b:mol.frame[i].getBonds())
				bonds.add(b);
				
		for (Bond b:bonds) 
			new Bond(getAtomNo(b.getA1().getNo()),getAtomNo(b.getA2().getNo()),b.getType());
		
		this.flipped = new HashSet<Bond>(mol.flipped);	
		
		refreshAtomNumbers(); // just in case ..
	}
	
	// for ring-closure OH-OH -> OLD
	public Molecule(Molecule mol,Atom O1,Atom H1,Atom O2,Atom H2) { // O2 and H1,H2 get lost as water
		//getting ready
		this.frame = new Atom[mol.frame.length-1];
		this.protons = new LinkedList<Atom>();
		this.name = new String(mol.name+"-cl");
		this.explHs = mol.explHs;
		
		int x = 0;
		//filling with NEW (!) objects
		for (int i=0;i<mol.frame.length;i++) {
			if (mol.frame[i]==O2){ // skipping O2
				x = 1;
				continue;
			}
			this.frame[i-x] = new Atom(mol.frame[i]);
			
			// the left oxygen is a little bit shifted
			if (this.frame[i-x].equals(O1)) {
				this.frame[i-x].setX( (O1.getX()-O2.getX())*0.5f+O2.getX() );
				this.frame[i-x].setY( (O1.getY()-O2.getY())*0.5f+O2.getY() );
				this.frame[i-x].setZ( (O1.getZ()-O2.getZ())*0.5f+O2.getZ() );
			}
		}
		
		for (Atom a:mol.getProtons()) {
			if (a==H1 || a==H2)
				continue;
			this.protons.add(new Atom(a));
		}
		
		//adding NEW (!) bonds
		Set<Bond> bonds = new HashSet<Bond>();
		for (int i=0;i<mol.frame.length;i++) {
			if (mol.frame[i]==O2)
				continue;
			for (Bond b:mol.frame[i].getBonds())
				bonds.add(b);
		}	
		for (Bond b:bonds) {
			if (b.getA1().equals(H1)||b.getA2().equals(H1)||b.getA1().equals(H2)||b.getA2().equals(H2))
				continue;
			
			else if (b.getA1().equals(O2) || b.getA2().equals(O2)) {
				new Bond(getEqualAtom(b.getOtherAtom(O2)),getEqualAtom(O1),b.getType());
			}
			else
				new Bond(getEqualAtom(b.getA1()),getEqualAtom(b.getA2()),b.getType());
			
		}
		this.flipped = new HashSet<Bond>(mol.flipped);	
		refreshAtomNumbers(); // Cs etc ...
	}
	
	// for ring-closure OH-O
		public Molecule(Molecule mol,Atom O1,Atom H1,Atom O2,boolean side) { // reconnect and replace!
			
			//getting ready
			this.frame = new Atom[mol.frame.length];
			this.protons = new LinkedList<Atom>();
			this.name = new String(mol.name+"-cl");
			this.explHs = mol.explHs;
			
			Atom C2 = O2.getAdjFrame(2,false).get(0);
			
			Atom O1n = null;
			Atom O2n = null;
			Atom H1n = null;
			Atom C2n = null;
			
			//filling with NEW (!) objects
			for (int i=0;i<mol.frame.length;i++) {
				this.frame[i] = new Atom(mol.frame[i]);
				
				if (this.frame[i].equals(O1))
					O1n = this.frame[i];
				else if (this.frame[i].equals(O2))
					O2n = this.frame[i];
				else if (this.frame[i].equals(C2))
					C2n = this.frame[i];
			}
			
			boolean skip = false;
			for (Atom a:mol.getProtons()) {
				if (!skip && a.equals(H1)) {
					H1n = new Atom(a);
					this.protons.add(H1n);
					skip = true;
				}
				else
					this.protons.add(new Atom(a));
			}
			
			
			//determine ALL old bonds
			Set<Bond> bonds = new HashSet<Bond>();
			for (int i=0;i<mol.frame.length;i++) 
				for (Bond b:mol.frame[i].getBonds())
					bonds.add(b);
			
			// rebuild new ones
			for (Bond b:bonds) {
				
				// this bond is set later!
				if (b.getA1().equals(H1)||b.getA2().equals(H1))
					continue;
				
				else if (b.getA1().equals(O2) || b.getA2().equals(O2)) {
					new Bond(C2n,O1n,(byte) 1); // new ring closure
				}
				else
					new Bond(getEqualAtom(b.getA1()),getEqualAtom(b.getA2()),b.getType());
				
			}
			
			// rough placement O
			O1n.setX( (O1.getX()-O2.getX())*0.5f+O2.getX() );
			O1n.setY( (O1.getY()-O2.getY())*0.5f+O2.getY() );
			O1n.setZ( (O1.getZ()-O2.getZ())*0.5f+O2.getZ() );
			// new bond is already set
			
			// placement =O
			positionH(O2n,C2n,side);
			new Bond(C2n,O2n,(byte) 1); // was double bond before
			// and adding the bond (afterwards to avoid distraction in placement and possibility to reuse already written function)
			
			// placement of H1 (formerly adjacent O1, not to O2(n))
			positionH(H1n,O2n);
			
			// add last bond
			new Bond(O2n,H1n,(byte) 1);
			
			this.flipped = new HashSet<Bond>();	
	}
	
	public Molecule(Molecule mol, int ring, int mid, int donor) {
		//getting ready
		this.frame = new Atom[mol.frame.length];
		this.protons = new LinkedList<Atom>();
		this.name = new String(mol.name);
		this.explHs = mol.explHs;

		//filling with NEW (!) objects
		for (int i=0;i<mol.frame.length;i++) 
			this.frame[i] = new Atom(mol.frame[i]);
		for (Atom a:mol.getProtons())
			this.protons.add(new Atom(a));

		//adding NEW (!) bonds
		Set<Bond> bonds = new HashSet<Bond>();
		for (int i=0;i<mol.frame.length;i++) 
			for (Bond b:mol.frame[i].getBonds())
				bonds.add(b);

		for (Bond b:bonds) 
			new Bond(getAtomNo(b.getA1().getNo()),getAtomNo(b.getA2().getNo()),b.getType());

		this.flipped = new HashSet<Bond>(mol.flipped);	
		
		// open-up
		Atom r = this.getAtomNo(ring);
		Atom m = this.getAtomNo(mid);
		Atom don = this.getAtomNo(donor);
		Atom H = don.getAdjHs().get(0);
		
		r.getBondTo(m).delBond();
		don.getBondTo(H).delBond();
		don.getBondTo(m).setType(2);
		
		new Bond(r, H, (byte) 1);
		positionH(H,r);
	}
	
	public Molecule(String file) {
		MolFileReader reader=new HinReader(file);
		
		// other possible formats ...
		
		frame = reader.getFrame();
		protons = reader.getProtons();
		name = reader.getName();
		explHs = reader.isExplHs();
		flipped = new HashSet<Bond>();
		refreshAtomNumbers();// important for array,getAtomNo,C-mapping!
	}
	
	public Molecule(InputStream instream) {
		MolFileReader reader=new StreamReader(instream);
		
		// other possible formats ...
		
		frame = reader.getFrame();
		protons = reader.getProtons();
		name = reader.getName();
		explHs = reader.isExplHs();
		flipped = new HashSet<Bond>();
		refreshAtomNumbers();
	}
	
	// for fragments
	public Molecule(Collection<Atom> atoms) {
		int n=0;
		
		for (Atom a:atoms)
			if (a.getType()!=Element.H)
				n++;
		
		this.protons = new LinkedList<Atom>();
		this.frame = new Atom[n];
		this.name = "fragment";
		Set<Bond> bonds = new HashSet<Bond>();
		
		n = 0;
		
		for (Atom a:atoms) {
			if (a.getType()==Element.H)
				protons.add(new Atom(a));
			else
				frame[n++] = new Atom(a);
			
			for (Bond b:a.getBonds())
				if (atoms.contains(b.getOtherAtom(a)))
					bonds.add(b);
		}
		
		for (Bond b:bonds) {
			Atom a1 = null;
			Atom a2 = null;
			
			for (Atom a:frame) {
				if (a.getNo()==b.getA1().getNo())
					a1 = a;
				else if (a.getNo()==b.getA2().getNo())
					a2 = a;
			}
			
			if (a1==null || a2==null)
				for (Atom a:protons) {
					if (a.getNo()==b.getA1().getNo())
						a1 = a;
					else if (a.getNo()==b.getA2().getNo())
						a2 = a;
				}
			
			// otherwise errors in SMILES ... adding Hs yields nice result
			if (a1 == null)
				new Bond(a2,new Atom(Element.H,n,0,0,0,0),(byte)1);
			else if (a2 == null)
				new Bond(a1,new Atom(Element.H,n,0,0,0,0),(byte)1);
			else
				new Bond(a1,a2,b.getType());
		}
		
		refreshAtomNumbers();
	}
	
	public static List<Molecule> getMoleculesFromFolder(String path) {
		List<Molecule> molecules = new LinkedList<Molecule>();
		
		FilenameFilter ff = new FilenameFilter() {
			
			@Override
			public boolean accept(File dir, String name) {
				
				return name.endsWith(".hin")||name.endsWith(".HIN");
			}
		};
		
		File[] files = new File(path).listFiles(ff);
		Arrays.sort(files);
		if (files == null)
			return molecules;
		
		Molecule temp = null;
		for (File f:files){
			temp = new Molecule(f.getAbsolutePath());
			temp.setName(f.getAbsolutePath());
			molecules.add(temp);
		}
		
		temp = null;
		
		return molecules;
	}
	
	public static void writeMoleculesToFolder(Collection<Molecule> molecules,String path) {
		int n = 0;
		if (!path.endsWith("/"))
			path += "/";
		for (Molecule m:molecules)
			m.writeHinFile(path+(n++)+".hin");
	}
	
	
	public List<Atom> getHDon() {
		List<Atom> donors = new LinkedList<Atom>();
		
		for (Atom a:frame)
			if (a.isHDonor())
				donors.add(a);
		
		return donors;
	}
	
	public List<Atom> getHAcc(boolean sp2_flag) {
		List<Atom> acceptors = new LinkedList<Atom>();
		 
		for (Atom a:frame)
			if (a.isHAcc(sp2_flag))
				acceptors.add(a);
		
		return acceptors;
	}
	
	public Set<Bond> getBonds() {
		Set<Bond> all_bonds = new HashSet<Bond>();
		
		for (Atom a:frame)
			all_bonds.addAll(a.getBonds());
		
		return all_bonds;
	}
	
	public void refreshAtomNumbers() {
		//Cs first
		for (int i=0;i<frame.length-1;i++) {
			if (frame[i].getType()==Element.C)
				continue;
			for (int j=i+1;j<frame.length;j++) {
				// swap
				if (frame[j].getType()==Element.C) {
					Atom temp = frame[i];
					frame[i] = frame[j];
					frame[j] = temp;
					break;
				}
			}
		}
		int n=1;
		for (Atom a:frame)
			a.setNo(n++);
		for (Atom a:protons)
			a.setNo(n++);
	}
	
	public void unMarkAtoms() {
		for (Atom a:frame)
			a.unMark();
		for (Atom a:protons)
			a.unMark();
		from = null;
	}
	
	public Atom getAtomNo(int i) {
		if (i-1<frame.length)
			return frame[i-1];
		else
			return protons.get(i-1-frame.length);
	}
	
	public Atom getEqualAtom(Atom a) {
		if (a.getType() != Element.H) {
			for (Atom f:frame)
				if (f.equals(a))
					return f;
		}
		else {
			for (Atom h:protons)
				if (h.equals(a))
					return h;
		}
		
		return null;
	}
	
	public void deleteH(Atom H) {
		for (Bond b:H.getBonds())
			b.delBond();
		protons.remove(H);
		refreshAtomNumbers();
	}
	
	public void addNewH(Atom at) {
		Atom H = new Atom(Element.H,Integer.MAX_VALUE,0,0,0,0);
		protons.add(H);
		positionH(H,at,false);
		new Bond(at,H,(byte)1);
		refreshAtomNumbers();
	}
	
	public String toString() {
		StringBuffer molecule = new StringBuffer();
		molecule.append("mol 1 \""+name+"\" \n");
		for (Atom a:frame) 
			molecule.append(a.toString()+"\n");
		for (Atom a:protons)
			molecule.append(a.toString()+"\n");
		
		molecule.append("endmol 1\n");
	
		return molecule.toString();
	}
	
	public void neutralizeCharges() {
		for (Atom a:frame)
			a.setCharge(0);
		for (Atom a:protons)
			a.setCharge(0);
	}
	
	/**
	 * Produces a map mapping H-donors to all its "reachable" H-acceptors
	 * @return Map (don)->(reachable acc)
	 */
	public Map<Atom,Set<Atom>> getFilteredAltReachMap(boolean aromatic_breakage, boolean sp2_flag) {
		Map<Atom,Set<Atom>> all = getAltReachMap(aromatic_breakage,false);
		Set<Atom> acceptors = new HashSet<Atom>(getHAcc(sp2_flag));
		Set<Atom> temp;
		
		//filtering for acceptors in the reachability lists
		Map<Atom,Set<Atom>> map=new HashMap<Atom, Set<Atom>>();
		
		for (Atom a:all.keySet()) {
			temp = new HashSet<Atom>();
			for (Atom b:all.get(a))
				if (acceptors.contains(b))
					temp.add(b);
			if (temp.size()>0)
				map.put(a,temp);
			
		}
		
		return map;
	}
	
	/**
	 * Produces a map mapping H-donors to all alternating reachable atoms
	 * @param aromatic_breakage
	 * @return
	 */
	public Map<Atom,Set<Atom>> getAltReachMap(boolean aromatic_breakage,boolean whole_path) {
		Map<Atom,Set<Atom>> all = new HashMap<Atom, Set<Atom>>();
		Deque<Atom> stack1 = new LinkedList<Atom>();
		Deque<Atom> stack2 = new LinkedList<Atom>();
		Set<Atom> temp, m1, m2;
		
		m1 = new HashSet<Atom>();
		m2 = new HashSet<Atom>();
		
		unMarkAtoms();
		List<Atom> donors = getHDon();
		
		
		for (Atom a:donors) {
			temp = new HashSet<Atom>();
			stack1.push(a);
			m1.clear();
			m2.clear();
			Atom current;
			while (stack1.size()>0 || stack2.size()>0) {
				
				while (stack1.size()>0){
					current = stack1.pop();
					m2.add(current);
					for (Atom b:current.getAdjFrame(1,aromatic_breakage))
						if (!m1.contains(b)) {
							stack2.add(b);
							if (whole_path)
								temp.add(b);
						}
				}
				
				while (stack2.size()>0) {
					current = stack2.pop();
					m1.add(current);
					for (Atom b:current.getAdjFrame(2,aromatic_breakage))
						if (!m2.contains(b)) {
							stack1.add(b);
							temp.add(b);
						}
				}
				
			}
			
			unMarkAtoms();
			all.put(a,temp);
		}
		
		return all;
	}
	
	public boolean isExplHs() {
		return explHs;
	}

	public void setExplHs(boolean explHs) {
		this.explHs = explHs;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Atom[] getFrame() {
		return frame;
	}

	public List<Atom> getProtons() {
		return protons;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (explHs ? 1231 : 1237);
		result = prime * result + Arrays.hashCode(frame);
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + ((protons == null) ? 0 : protons.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) { // TODO molecule hashcode/equals
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof Molecule))
			return false;
		Molecule other = (Molecule) obj;
		if (explHs != other.explHs)
			return false;
		if (!Arrays.equals(frame, other.frame))
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (protons == null) {
			if (other.protons != null)
				return false;
		} else if (!protons.equals(other.protons))
			return false;
		
		try {
			for (int i=0;i<frame.length;i++)
				if (!frame[i].equals(other.frame[i]))
					return false;
		} catch (Exception E) {
			return false;
		}
		return true;
	}
	
	public void transform(int donor,int acceptor,List<Integer> bond_list) {
		
		// donor
		Atom don = getAtomNo(donor);
		don.resetAtomType();
		
		// acceptor has exactly one db/tb
		Atom acc = getAtomNo(acceptor);
		acc.resetAtomType();
		
		// change bonds between
		for (Bond b:getBondsBetween(bond_list)) {
				
			if (!flipped.contains(b)) {
				b.flip();
				flipped.add(b);
			}
		}
		
		if (explHs) {
			// choose the H
			Atom H = don.getAdjHs().get(0);
			// reset other Hs?
		     
			// move the H
			positionH(H,acc);
			H.resetAtomType();
			// bind it
			don.getBondTo(H).setBindToInsteadOf(acc,don);
		}
		
	}
	
	/**
	 * Places an unconnected hydrogen (or other atom :-) )
	 * @param H to place
	 * @param acc H-acceptor
	 */
	public void positionH(Atom H,Atom acc) {
		positionH(H, acc, false);
	}
	
	/**
	 * Places an unconnected hydrogen (or other atom :-) )
	 * @param H to place
	 * @param acc H-acceptor
	 * @param reverse allows for usage in stereo-problems
	 */
	public void positionH(Atom H,Atom acc,boolean reverse) {
		List<Atom> adj = acc.getAdjAtoms();
		Atom next1,next2,next3;
		float x,y,z,factor;
		switch (adj.size()) {
		case 1:
			next1 = adj.get(0);
			H.setX(1.8f*(acc.getX()-next1.getX())+next1.getX()); 
			H.setY(1.8f*(acc.getY()-next1.getY())+next1.getY()); 
			H.setZ(1.8f*(acc.getZ()-next1.getZ())+next1.getZ());
			return;
		case 2:
			next1 = adj.get(0);
			next2 = adj.get(1);
			x = next1.getX()-acc.getX() + next2.getX()-acc.getX();
			y = next1.getY()-acc.getY() + next2.getY()-acc.getY();
			z = next1.getZ()-acc.getZ() + next2.getZ()-acc.getZ();
			factor = -1.0f/(float) Math.sqrt(x*x+y*y+z*z);
			x *= factor;
			y *= factor;
			z *= factor;
			H.setX(acc.getX()+x);
			H.setY(acc.getY()+y);
			H.setZ(acc.getZ()+z);
			return;
		default: // case 3+
			next1 = adj.get(0);
			next2 = adj.get(1);
			next3 = adj.get(2);
			
			// u2*v3 - u3*v2
			x = (next1.getY()-next3.getY())*(next2.getZ()-next3.getZ()) - (next1.getZ()-next3.getZ())*(next2.getY()-next3.getY());
				
			// u3*v1 - u1*v3
			y = (next1.getZ()-next3.getZ())*(next1.getX()-next3.getX()) - (next1.getX()-next3.getX())*(next2.getZ()-next3.getZ());
			
			// u1*v2 - u2*v1
			z = (next1.getX()-next3.getX())*(next2.getY()-next3.getY()) - (next1.getY()-next3.getY())*(next2.getX()-next3.getX());
			
			// on which side of the hyperplane is acc
			float dot = acc.getX()*x+acc.getY()+y+acc.getZ()*z;
			factor = (float)(Math.signum(dot)/Math.sqrt(x*x+y*y+z*z));
			
			// switches side			
			if (reverse)
				factor *= -1;
			
			x *= factor;
			y *= factor;
			z *= factor;
			H.setX(acc.getX()+x);
			H.setY(acc.getY()+y);
			H.setZ(acc.getZ()+z);
			return;
			
		}
	}
	
	/**
	 * Replaces (for stereo-conv) an connected hydrogen
	 * @param H to place
	 * @param acc H-acceptor
	 */
	public void inverseH(Atom H,Atom acc) {
		float x = acc.getX() - H.getX();
		float y = acc.getY() - H.getY();
		float z = acc.getZ() - H.getZ();
		H.setX(acc.getX()+x);
		H.setY(acc.getY()+y);
		H.setZ(acc.getZ()+z);
	}
	
	public List<Molecule> getOtherStereoForms() {
		List<Molecule> forms = new LinkedList<Molecule>();
		Molecule temp;
		
		for (Atom a:getChiralAtoms()) {
			if (a.getAdjHs().size() == 1) {
				int atom_no = a.getAdjHs().get(0).getNo();
				temp = new Molecule(this);
				Atom invertible = temp.getAtomNo(atom_no);
				temp.inverseH(invertible, invertible.getAdjAtoms().get(0));
				forms.add(temp);
			}
		}
		
		return forms;
	}
	
	public void writeHinFile(String file) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			writer.write("mol 1 \""+name+"\"");
			writer.newLine();
			
			for (Atom a:frame) {
				writer.write(a.toString());
				writer.newLine();
				}
			for (Atom a:protons) {
				writer.write(a.toString());
				writer.newLine();
				}
			
			writer.write("endmol 1");
			writer.newLine();
			writer.flush();
			writer.close();
		} catch (IOException e) {
			System.out.println("Problem writing "+file+".");
			System.exit(1);
		}
	}
	
	public List<Atom> getChiralAtoms() {
		
		if (chiral_atoms == null)
			chiral_atoms = BabelInterface.get_chiral_atoms(this);
		
		return chiral_atoms;
	}
	
	public Bond[] getBondsBetween(List<Integer> list) {
		Bond[] bonds = new Bond[list.size()-1];
		Integer[] atoms = list.toArray(new Integer[list.size()]);
		
		for (int i=0;i<atoms.length-1;i++)
			bonds[i] = (getAtomNo(atoms[i]).getBondTo(getAtomNo(atoms[i+1])));
		
		return bonds;
	}
	
	public boolean isPlausible() {
		for (Atom a:frame)
			if (!a.isPlausible()) {
				//System.out.println(a);
				return false;
			}
		return true;
	}
	
	public void detectCycleAtoms() {
		if (cycle_atoms != null)
			return;
		
		from = new HashMap<Atom, Atom>();
		cycle_atoms = new HashSet<Atom>();
		
		// dfs + cycle-detection -> filling of cycle_bonds 
		
		for (Atom current:frame) {
			if (current.isMarked())
				continue;
			cDFS(current);
		}
		
		unMarkAtoms();
	}
	
	public Set<Atom> getCycleAtoms() {
		
		detectCycleAtoms();
		
		return cycle_atoms;
	}
	
	public Map<Bond,List<Atom>> getRotatableBondsMap() {
		
		detectCycleAtoms();
		
		List<Bond> rotatable = new LinkedList<Bond>();
		
		for (Bond b:getBonds()) {
			
			// bonds with Hs
			if (b.getA1().getType()==Element.H || b.getA2().getType()==Element.H)
				continue;
			
			// non-singlebonds
			if (b.getType()!=1)
				continue;
			
			// bonds in cycles
			if (cycle_atoms.contains(b.getA1()) && cycle_atoms.contains(b.getA2()))
				continue;
			
			// more like partial db
			if ((b.getA1().getType()==Element.N && b.getA2().getType()==Element.C) ||  // with this rule
					(b.getA1().getType()==Element.C && b.getA2().getType()==Element.N)) // all resonance forms of peptide bonds are hit
				continue;
			
			rotatable.add(b);
		}
		
		from = null;
		cycle_atoms = null;
		
		Map<Bond,List<Atom>> rot_map = new HashMap<Bond, List<Atom>>();
		
		int sum;
		List<Atom> temp;
		for (Bond b:rotatable) {
			
			temp = b.getMinimalRotAtoms();
			sum = 0;
			for (Atom a:temp) {
				if (a.getType()==Element.H)
					continue;
				sum++;
			}
			
			// stuff like methyl doesn't need to be rotated
			if (sum>1)
				rot_map.put(b,temp);
			unMarkAtoms();
		}
		
		return rot_map;
	}
	
	private void cDFS(Atom a) {
		a.mark();
		for (Atom adj:a.getAdjFrame(0,true)) {
			if (!adj.isMarked()) {
				from.put(adj,a);
				cDFS(adj);
			}
			else if (from.get(a)!= null && !from.get(a).equals(adj)) { // not the edge returning
				// cycle-backtrack
				if (!cycle_atoms.contains(a)){
					Atom temp = a;
					while (temp!=adj && temp!=null) {
						cycle_atoms.add(temp);
						temp = from.get(temp);
					}
				cycle_atoms.add(temp);
				}
			}
		}
	}
	
	public static Map<Bond,Set<Atom>> getDisjunctRotAtomsMap(Map<Bond,List<Atom>> rot_map) {
		Map<Bond,Set<Atom>> disj = new HashMap<Bond, Set<Atom>>();
		
		// generate sets
		for (Bond b:rot_map.keySet())
			disj.put(b,new HashSet<Atom>(rot_map.get(b)));
		
		Set<Atom> set_a;
		for (Bond a:rot_map.keySet()){
			set_a = disj.get(a);
			for (Bond b:rot_map.keySet()) {
				if (a==b || rot_map.get(a).size()<rot_map.get(b).size())
					continue;
				set_a.removeAll(rot_map.get(b));
			}
		}
		
		return disj;
	}
	
	public List<Atom> getFixAtoms(Map<Bond,Set<Atom>> disj_rot_atoms) {
		Set<Atom> fix = new HashSet<Atom>(protons);
		for (Atom a:frame)
			fix.add(a);
		
		for (Set<Atom> s:disj_rot_atoms.values())
			fix.removeAll(s);
		
		return new LinkedList<Atom>(fix);
	}
	
	public List<Atom> getAllAtoms() {
		List<Atom> fix = new LinkedList<Atom>();
		
		for (Atom a:frame)
			fix.add(a);
		
		fix.addAll(protons);
		return fix;
	}
	
	public String getInstance() {
		return super.toString();
	}
	
	public void buildFastScoringMap() {
		score_map = new HashMap<Atom, List<Atom>>();
		score_map_4 = new HashMap<Atom, List<Atom>>();
		
		Atom a;
		Set<Atom> set;
		List<Atom> temp;
		List<Atom> temp_4;
		for (int i=1;i<=frame.length+protons.size();i++) {
			a = getAtomNo(i);
			set = new HashSet<Atom>(getAllAtoms());
			temp = new LinkedList<Atom>();
			temp_4 = new LinkedList<Atom>();
			set.remove(a);
			for (Atom b:a.getAdjAtoms()) {
				set.remove(b);
				for (Atom c:b.getAdjAtoms()) {
					if (c==a)
						continue;
					set.remove(c);
					for (Atom d:c.getAdjAtoms()) {
						if (d==b || d==a)
							continue;
						if (a.getNo()<d.getNo())
							temp_4.add(d);
						set.remove(d);
					}
						
				}
			}
			
			for (Atom b:set) {
				if (a.getNo()<b.getNo())
					temp.add(b);
			}
			set = null;
			score_map.put(a,temp);
			score_map_4.put(a,temp_4);
		}
	}
	
	public double getFastScore(double E_0) {
		double score = 0;
		double r;
		if (E_0<=0)
			E_0 = 1.0;
		
		if (score_map==null)
			buildFastScoringMap();
		
		Atom a;
		double t1,t2;
		for (int i=1;i<=frame.length+protons.size();i++) {
			a = getAtomNo(i);
			for (Atom b:score_map.get(a)) {
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
				score += 331.92552431411053 * a.getCharge()*b.getCharge()/ 
				( r*E_0);	
			}
			
			// 1-4 only *0.5
			for (Atom b:score_map_4.get(a)) {
				r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
						+(a.getY()-b.getY())*(a.getY()-b.getY())
						+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
				
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
		return score;
	}
	
	public double getCoulomb(double E_0) {
		double score = 0.0;
		double r;
		if (E_0<=0)
			E_0 = 1.0;
		
		if (score_map==null)
			buildFastScoringMap();
		
		Atom a;
		for (int i=1;i<=frame.length+protons.size();i++) {
			a = getAtomNo(i);
			for (Atom b:score_map.get(a)) {
				r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
						+(a.getY()-b.getY())*(a.getY()-b.getY())
						+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
				
				r += 0.05;
				
				// coulomb in kcal / mol
				score += a.getCharge()*b.getCharge()/ r;	
			}
			
			// 1-4 only *0.75
			for (Atom b:score_map_4.get(a)) {
				r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
						+(a.getY()-b.getY())*(a.getY()-b.getY())
						+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
				
				r += 0.05;
				
				// coulomb in kcal / mol
				score += 0.75 * a.getCharge()*b.getCharge() / r;	
			}
		}
		return score*332.0716/E_0;
	}
	
	public double getVdW() {
		double score = 0;
		double r;
		
		if (score_map==null)
			buildFastScoringMap();
		
		Atom a;
		double t1,t2;
		for (int i=1;i<=frame.length+protons.size();i++) {
			a = getAtomNo(i);
			for (Atom b:score_map.get(a)) {
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
			}
			
			// 1-4 only *0.5
			for (Atom b:score_map_4.get(a)) {
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
			}
		}
		return score;
	}
	
	public double getsoftVdW() {
		double score = 0;
		double r;
		
		if (score_map==null)
			buildFastScoringMap();
		
		Atom a;
		double t1,t2;
		for (int i=1;i<=frame.length+protons.size();i++) {
			a = getAtomNo(i);
			for (Atom b:score_map.get(a)) {
				r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
						+(a.getY()-b.getY())*(a.getY()-b.getY())
						+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
				
				r += 0.05;
				// VdW
				t1 = (a.getVdW()+ b.getVdW())/r;
				double term6 = t1 * t1 * t1;
				term6 = term6 * term6;
				double term12 = term6 * term6;
				
				t2 = Math.sqrt(a.getVdWEpsilon()*b.getVdWEpsilon());
				score += t2 * (term12 - 2*term6);
			}
			
			// 1-4 only *0.5
			for (Atom b:score_map_4.get(a)) {
				r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
						+(a.getY()-b.getY())*(a.getY()-b.getY())
						+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
				
				r += 0.05;
				// VdW
				t1 = (a.getVdW()+ b.getVdW())/r;
				double term6 = t1 * t1 * t1;
				term6 = term6 * term6;
				double term12 = term6 * term6;
				
				t2 = Math.sqrt(a.getVdWEpsilon()*b.getVdWEpsilon());
				score += 0.5 * t2 * (term12 - 2*term6);
			}
		}
		return score;
	}
	
	public double getTruncVdW() {
		double score = 0;
		double r;
		
		if (score_map==null)
			buildFastScoringMap();
		
		Atom a;
		double t1,t2;
		for (int i=1;i<=frame.length+protons.size();i++) {
			a = getAtomNo(i);
			for (Atom b:score_map.get(a)) {
				r = Math.sqrt((a.getX()-b.getX())*(a.getX()-b.getX())
						+(a.getY()-b.getY())*(a.getY()-b.getY())
						+(a.getZ()-b.getZ())*(a.getZ()-b.getZ()));
				
				// VdW
				t1 = (a.getVdW()+ b.getVdW())/r;
				double term6 = t1 * t1 * t1;
				term6 = term6 * term6;
				double term12 = term6 * term6;
				
				t2 = Math.sqrt(a.getVdWEpsilon()*b.getVdWEpsilon());
				
				// truncation
				if (t2 * (term12 - 2*term6)>0)
					score += t2 * (term12 - 2*term6);
			}
			
			// 1-4 only *0.5
			for (Atom b:score_map_4.get(a)) {
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
			}
		}
		return score;
	}
	
	public static double getTorsionAngle(Atom a1,Atom a2,Atom a3,Atom a4) {
		float[] vec1 = new float[3];
		float[] vec2 = new float[3];
		float x1,x2,x3,y1,y2,y3;
	
		// normalvector of first plane
		x1 = a1.getX() - a2.getX(); x2 = a1.getY() - a2.getY(); x3 = a1.getZ() - a2.getZ();
		y1 = a3.getX() - a2.getX(); y2 = a3.getY() - a2.getY(); y3 = a3.getZ() - a2.getZ();
		vec1[0] = x2*y3 - x3*y2;
		vec1[1] = x3*y1 - x1*y3;
		vec1[2] = x1*y2 - x2*y1;

		// normalvector of second plane
		x1 = a2.getX() - a3.getX(); x2 = a2.getY() - a3.getY(); x3 = a2.getZ() - a3.getZ();
		y1 = a4.getX() - a3.getX(); y2 = a4.getY() - a3.getY(); y3 = a4.getZ() - a3.getZ();
		vec2[0] = x2*y3 - x3*y2;
		vec2[1] = x3*y1 - x1*y3;
		vec2[2] = x1*y2 - x2*y1;

		// calculate angle between these vectors
		double angle = Math.toDegrees(Math.acos((vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]) / 
				(Math.sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])
						*Math.sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]))));
			
		return angle;
	}
	
	public List<Set<Atom>> getFragments() {
		List<Set<Atom>> fragments = new LinkedList<Set<Atom>>();
		Set<Bond> rot_bond = getRotatableBondsMap().keySet();
		Set<Atom> rot_bond_atoms = new HashSet<Atom>();
		
		// store all atoms belonging to rotatable bonds
		for (Bond b:rot_bond) {
			rot_bond_atoms.add(b.getA1());
			rot_bond_atoms.add(b.getA2());
		}
		
		unMarkAtoms();
		Set<Atom> current_fragment = null;
		Deque<Atom> stack = new LinkedList<Atom>();
		Atom current = null;
		for (Atom start:rot_bond_atoms) {
			
			// every atom at a rotatable bond can only belong to one fragment
			if (start.isMarked())
				continue;
			
			current_fragment = new HashSet<Atom>();
			stack.push(start);
			
			while (stack.size()>0) {
				current = stack.pop();
				if (current.isMarked())
					continue;
				current.mark();
			
				current_fragment.add(current);
				
				for (Bond b:current.getBonds()) {
					if (rot_bond.contains(b))
						continue;
					stack.push(b.getOtherAtom(current));
				}
			}
			
			fragments.add(current_fragment);
		}
		
		unMarkAtoms();
		
		return fragments;
	}

	public void CheckAromaticSixRingCs() {
		Queue<Atom> SixRingCs = new LinkedList<Atom>();
		for (Atom a:frame) {
			if (a.getType()!=Element.C)
				break;
			if (a.isInAromaticSixMemberRing())
				SixRingCs.add(a);
		}
		
		if (SixRingCs.size() < 6)
			return;
		
		// order ones with double bonds first (bordercases!)
		Atom[] temp = new Atom[SixRingCs.size()];
		int i = 0;
		while (SixRingCs.size()>0) {
			temp[i++] = SixRingCs.remove();
		}
		
		for (i = 0;i<temp.length-1;i++) {
			for (int j = i+1;j<temp.length;j++) {
				int[] sum_i=new int[4];
				// count different types of bonds
				for (Bond b:temp[i].getBonds()) 
					sum_i[b.getType()-1]++;
			
				int[] sum_j=new int[4];
				// count different types of bonds
				for (Bond b:temp[j].getBonds()) 
					sum_j[b.getType()-1]++;
				
				if (sum_j[1]==1 && sum_i[1]==0) { // wer doppelbindungen hat kommt vor!
					Atom bla = temp[i];
					temp[i] = temp[j];
					temp[j] = bla;
					break; // temp[i] besetzt...
				}
			}
		}
		
		int db = 0;
		for (Bond b:temp[0].getBonds())
			if (b.getType() == 2)
				db++;
		
		if (db == 0) // nothing to fix!
			return;
				
		for (i = 0;i<temp.length;i++) {
			SixRingCs.add(temp[i]);
		}
		
		// check
		Atom current;
		while (SixRingCs.size()>0) {
			current = SixRingCs.remove();
			
			int[] sum=new int[4];
			// count different types of bonds
			for (Bond b:current.getBonds()) 
				sum[b.getType()-1]++;
					
			// bordercase
			if (sum[1]==1 && sum[3]==2) {
				for (Bond b:current.getBonds()) 
					if (b.getType()==4)
						b.setType(1); // they become single bonds
				continue;
			}
			
			// middle case
			if (sum[1]==0 && sum[3]==1 && sum[0]==2) {
				for (Bond b:current.getBonds()) 
					if (b.getType()==4)
						b.setType(2); // they become double bonds
				continue;
			}
			
			// middle case
			if (sum[1]==0 && sum[3]==1 && sum[0]==3) {
				for (Bond b:current.getBonds()) 
					if (b.getType()==4)
						b.setType(1); // they become single bonds
				continue;
			}
				
			
			// other case
			if (sum[1]==0 && sum[3]==2 && sum[0]==2) {
				for (Bond b:current.getBonds()) 
					if (b.getType()==4)
						b.setType(1); // they become single bonds
				continue;
			}
			
			// other case 2
			if (sum[1]==0 && sum[3]==2 && sum[0]==1) {
				Bond[] bonds = new Bond[2];
				i = 0;
				for (Bond b:current.getBonds()) 
					if (b.getType()==4)
						bonds[i++] = b;
				
				for (Bond b:bonds) {
					Atom other = b.getOtherAtom(current);
					int[] other_sum=new int[4];
					// count different types of bonds
					for (Bond bo:other.getBonds()) 
						other_sum[bo.getType()-1]++;
					if (other_sum[1]==0 && other_sum[3]==2 && other_sum[0]==1)
						b.setType(2); // the aromatic bond between them has to be aromatic!
				}
				continue;
			}
			
		}
		
	}
	
	public void straightenOH() {
		for (Atom a:frame)
			if (a.getType()==Element.O && a.getAdjHs().size()==1) {
				// hydrogen placement
				Atom H = a.getAdjHs().get(0);
				Atom next1 = a.getAdjFrame(0,false).get(0);
				H.setX(1.6f*(a.getX()-next1.getX())+next1.getX()); 
				H.setY(1.6f*(a.getY()-next1.getY())+next1.getY()); 
				H.setZ(1.6f*(a.getZ()-next1.getZ())+next1.getZ());
			}
	}
	
	public void resetAromatic() {
		for (Bond b:getBonds())
			if (b.getType() == 4)
				b.setType(1);
	}
	
	public boolean QaDValencyFix(Bond set) {
		
		// set bond set so the clauses are fulfilled
		if (set != null) {
			set.setType(2);
		}
		
		// formulate as CNF and solve with dpll
		List<List<Bond>> to_set = new LinkedList<List<Bond>>();
		
		for (Atom a:frame) {
			Float f = a.getValencyDiff();
			if (f == -1.0f) {
				List<Bond> temp = new LinkedList<Bond>();
				for (Bond b:a.getBonds()) {
					if (b.getOtherAtom(a).getValencyDiff() == -1.0f)
						temp.add(b);
				}
				to_set.add(temp);
			}
			if (f == 1.0f)
				return false; // probably later ...
		}
		
		// dpll
		// success: no poblems anymore
		if (to_set.size() == 0)
			return true;
		
		// fail-test and UP
		List<List<Bond>> to_set_temp = new LinkedList<List<Bond>>(to_set);
		List<Bond> to_del = new LinkedList<Bond>();
		
		for (List<Bond> l:to_set) {
			if (l.size() == 0) // empty
				return false;
			if (l.size() == 1)
				to_del.add(l.get(0));
		}
		
		// set bond and delete all clauses containing it
		for (Bond b:to_del) {
			b.setType(2);
			for (List<Bond> l:to_set)
				if (l.contains(b))
					to_set_temp.remove(l);
		}
		to_set = to_set_temp;
		
		// success: no poblems anymore
		if (to_set.size() == 0)
			return true;
		
		// choice
		for (Bond b:to_set.get(0)) {
			if (!QaDValencyFix(b))
				b.setType(1);
			else
				return true;
		}
		
		return false;
	}
	
	public boolean QaDValencyFix() {
		return QaDValencyFix(null);
	}
	
	@Override
	public int compareTo(Molecule o) {
		String[] a1 = name.split("/");
		String[] a2 = o.name.split("/");
		
		if (a1.length>1 && a2.length>1) {
			String n1 = a1[a1.length-1];
			String n2 = a2[a2.length-1];
			int x1 = Integer.parseInt(n1.split("-")[0]);
			int x2 = Integer.parseInt(n2.split("-")[0]);
			
			if (x1>x2)
				return +1;
			if (x1==x2)
				return 0;
			return -1;
		}
		
		return name.compareTo(o.name);
	}
	
}
