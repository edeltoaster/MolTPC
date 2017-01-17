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
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class Atom {
	private final Element type;
	private float x,y,z;
	private float charge;
	private int no;
	private boolean marked;
	private Set<Bond> bonds = new HashSet<Bond>();
	private float vdw = Float.NaN;
	private float vdw_eps = Float.NaN;
	private boolean isOH = false;
	private boolean isOHset = false;
	private boolean isKeto = false;
	private boolean isKetoSet = false;
	
	public Atom(Atom a) {
		this.type = a.type;
		this.no = a.no;
		this.x = a.x;
		this.y = a.y;
		this.z = a.z;
		this.charge = a.charge;
		//bonds done externally
	}
	
	public Atom(Element type,int no,float x,float y,float z,float charge) {
		this.type = type;
		this.no = no;
		this.x = x;
		this.y = y;
		this.z = z;
		this.charge = charge;
		//bonds are done externally
	}
	
	public void addBond(Bond b) {
		bonds.add(b);
	}
	
	public void deleteBond(Bond b) {
		bonds.remove(b);
	}
	
	public Set<Bond> getBonds() {
		return bonds;
	}
	
	/**
	 * Returns the matching bond
	 * @param a atom to which this bond binds to
	 * @return the bond or null if it does not exit
	 */
	public Bond getBondTo(Atom a) {
		for (Bond b:bonds)
			if (b.getOtherAtom(this).equals(a))
				return b;
		
		return null;
	}
	/**
	 * 
	 * @return could be donor?
	 */
	public boolean isHDonor() {
		int[] sum=new int[4];
		boolean H_seen = false;
		
		// count different types of bonds and adjacent Hs
		for (Bond b:bonds) {
			if (b.getType()==1)
				if (b.getOtherAtom(this).getType()==Element.H)
					H_seen = true;
			sum[b.getType()-1]++;
		}
		
		if (!H_seen)
			return false;
		
		if (sum[3]>0) {
			// sum all shared electrons, aromatic bonds counted as 1
			return (type.getValency() == sum[0]+sum[3]);
		}
		else if (sum[0] == bonds.size()) 
			return true; // restrict carbons connected to heteroatoms?
		
		return false;
	}
	
	/**
	 * 
	 * @return
	 */
	public boolean isHAcc(boolean sp2_flag) {
		int[] sum=new int[4];
		
		// count different types of bonds
		for (Bond b:bonds)
			sum[b.getType()-1]++;
		
		// sum all shared electrons, aromatic bonds counted as 1.5
		float bondsum = sum[0]+sum[1]*2+sum[2]*3+sum[3]*1.5f;
		
		// special dealing if there are aromatic bonds
		if (sum[3]>0) { // sp2_flag = true -> preservation
			if (sp2_flag && type==Element.C && isAromatic()) {// check: is in arom 6 ring? see tauthor!
				return false;
			}
			
			// now, explicit Hs given, one bond must have db character
			if ((type.getValency() > bondsum) && type==Element.C) {
				Atom d_adj = getAdjFrame(2,false).get(0);
				if (d_adj.getType() != Element.C && d_adj.getType() != Element.H)
					return false;
			} 
//			else if ((type.getValency() == bondsum) && type==Element.C) {
//				for (Atom a:getAdjFrame(4,false))
//					if (a.getType() != Element.C && a.getType() != Element.H)
//						return false;
//			}
			
			return (type.getValency() == bondsum);
			//return false;
		}
		// acceptors have exactly one db/tb
		else if ((sum[1]+sum[2]) == 1){
			//filter things like C=O -> C wont be acceptor
			if (type==Element.C) {
				Atom d_adj = getAdjFrame(2,false).get(0);
				if (d_adj.getType() != Element.C && d_adj.getType() != Element.H)
					return false;
			}
			
			return true;
		}
		// others don't take protons
		return false;
	}
	
	public List<Atom> getAdjHs() {
		List<Atom> adjHs=new LinkedList<Atom>();
		
		for (Bond b:bonds) {
			Atom temp = b.getOtherAtom(this);
			if (temp.getType() == Element.H)
				adjHs.add(temp);
		}
		
		return adjHs;
	}
	
	/**
	 * gets adjacent framework-atoms that are
	 * mode	0: connected
	 * 		1: connected with a singlebond
	 * 		2: connected with double/triple-bond
	 * @param mode mode of connection
	 * @return List of atoms s.a.
	 */
	public List<Atom> getAdjFrame(int mode,boolean aromatic_breakage) {
		List<Atom> frame = new LinkedList<Atom>();
		
		for (Bond b:bonds) {
			if (mode==0 || 
					(mode==1 && (b.getType()==1 || (b.getType()==4 && aromatic_breakage))) || 
					(mode==2 && (b.getType()==2 || b.getType()==3 || (b.getType()==4 && aromatic_breakage)))
					|| (mode==4 && b.getType()==4)) {
				if (b.getOtherAtom(this).getType() == Element.H)
					continue;
				if (mode != 0 && b.getType() == 4) { // reject sp3-N-paths etc
					int[] sum=new int[4];
					
					boolean H_seen = false;
					// count different types of bonds
					for (Bond b1:bonds) {
						sum[b1.getType()-1]++;
						if (b1.getOtherAtom(this).getType() == Element.H)
							H_seen = true;
					}
					
					if (!H_seen) {
						if (type.getValency() == sum[0]+sum[3])
							continue;
					}
				}
				frame.add(b.getOtherAtom(this));
			}
		}
		
		return frame;
	}
	
	public List<Atom> getAdjAtoms() {
		List<Atom> all = getAdjFrame(0,true);
		all.addAll(getAdjHs());
		return all;
	}
	
	public boolean isAdjacent(Atom a) {
		for (Bond b:bonds)
			if (b.getOtherAtom(this).equals(a))
				return true;
		return false;
	}
	
	public void rotate(Quaternion rot_q,Quaternion rot_q_conj,float[] corr) {
		Quaternion p = new Quaternion(x-corr[0],
				y-corr[1],
				z-corr[2]);
	
		float[] new_p = rot_q.mult(p).point_mult(rot_q_conj);
	
		x = (new_p[0]+corr[0]);
		y = (new_p[1]+corr[1]);
		z = (new_p[2]+corr[2]);
	}
	
	public List<Atom> getRotAtoms(Bond bond) {
		List<Atom> atoms = new LinkedList<Atom>();
		Deque<Atom> stack = new LinkedList<Atom>();
		
		// don't go in this direction
		bond.getOtherAtom(this).mark();
		
		stack.add(this);
		
		// BFS
		Atom current,next;
		while (stack.size()>0) {
			current = stack.pop();
			for (Bond b:current.getBonds()) {
				next = b.getOtherAtom(current);
				if (next.isMarked())
					continue;
				next.mark();
				stack.push(next);
				atoms.add(next);
			}
		}
		
		return atoms;
	}
	
	public Element getType() {
		return type;
	}

	public float getX() {
		return x;
	}

	public void setX(float x) {
		this.x = x;
	}

	public float getY() {
		return y;
	}

	public void setY(float y) {
		this.y = y;
	}

	public float getZ() {
		return z;
	}

	public void setZ(float z) {
		this.z = z;
	}

	public float getCharge() {
		return charge;
	}

	public void setCharge(float charge) {
		this.charge = charge;
	}

	public int getNo() {
		return no;
	}

	public void setNo(int no) {
		this.no = no;
	}
	
	public String toString() {
		String atom = "atom "+no+" - "+type.name()+" ** - "+charge+" "+x+" "+y+" "+z+" "+bonds.size();
		
		for (Bond b:bonds)
			atom += " "+b.getOtherAtom(this).getNo()+" "+b.getStrType();
		
		return atom;
	}
	
	public boolean isMarked() {
		return marked;
	}
	
	public void mark() {
		marked = true;
	}
	
	public void unMark() {
		marked = false;
	}
	
	public boolean onlySB() {
		int[] sum=new int[4];
		
		// count different types of bonds 
		for (Bond b:bonds) 
			sum[b.getType()-1]++;
	
		if (sum[3]>0) {
			// sum all shared electrons, aromatic bonds counted as 1
			int bondsum = sum[0]+sum[1]*2+sum[2]*3+sum[3];
			// now, explicit Hs given, one bond must have db character
			return (type.getValency()==bondsum);
		}
		else if (sum[0]==bonds.size()) {
			return true;
		}
		
		return false;
	}
	
	public boolean isPlausible() {
		float sum = 0;
		boolean arom = false;
		for (Bond b:bonds)
			if (b.getType()==4) {
				sum += 1.5;// ok?
				arom = true;
			}
			else
				sum += b.getType();
		
		// TODO plausibility for charged atoms ?
		
		if (arom)
			return (sum >= type.getValency());
		
		return (sum == type.getValency());
	}
	
	public float getValencyDiff() {
		float sum = 0;
		
		for (Bond b:bonds)
			if (b.getType()==4) {
				sum += 1.5;// ok?
			}
			else
				sum += b.getType();
		
		return (sum-type.getValency());
	}
	
	@Override
	public int hashCode() {
		return no; //fast hashCode needed
	}
	
	public boolean isOH() {
		if (isOHset)
			return isOH;
		
		if (type==Element.O) {
			isOH = (getAdjHs().size()!=0);
			isOHset = true;
			return isOH;
		}
		else 
			isOHset = true;
		
			return false;
	}
	
	public boolean isKeto() {
		if (isKetoSet)
			return isKeto;
		
		if (type==Element.O) {
			isKeto = (getAdjFrame(2,false).size() == 1);
			isKetoSet = true;
			return isKeto;
		}
		else
			isKetoSet = true;
		
			return false;
	}
	
	public boolean isNH2X() {
		if (type!=Element.N)
			return false;
		if (getAdjHs().size()<2)
			return false;
		
		return true;
	}
	
	public boolean isCOOH() {
		if (type!=Element.O)
			return false;

		if (getAdjHs().size()!=1)
			return false;

		if (getAdjFrame(1,false).size()!=1)
			return false;

		Atom C = getAdjFrame(1,false).get(0);
	
		if (C.getAdjFrame(2,false).size()!=1)
			return false;
		
		if (C.getAdjFrame(2,false).get(0).getType()!=Element.O)
			return false;
		
		return true;
	}
	
	public float getVdW() {
		if (Float.isNaN(vdw))
			vdw = type.getVdW(this);
		return vdw;
	}
	
	public float getVdWEpsilon() {
		if (Float.isNaN(vdw_eps))
			vdw_eps = type.getEps(this);
		return vdw_eps;
	}
	
	public void resetAtomType() {
		vdw = Float.NaN;
		vdw_eps = Float.NaN;
		isOHset = false;
	}
	
	public boolean isElWithdraw() {
		if (type==Element.H)
			return false;
		else if (type==Element.C) {
			for (Atom a:getAdjAtoms())
				if (a.type!=Element.C && a.type!=Element.H)
					return true;
			return false;
		}
		else
			return true;
	}
	
	public boolean isAromatic() {
		for (Bond b:bonds)
			if (b.getType()==4)
				return true;
		return false;
	}
	
	public boolean isInAromaticSixMemberRing() {
		Deque<Atom> stack = new LinkedList<Atom>();
		Map<Atom,Atom> came_from = new HashMap<Atom, Atom>();
		
		int steps = 0;
		stack.add(this);
		
		Atom current;
		while (stack.size()>0) {
			current = stack.pop();
			if ((current==this && steps>0) || steps>8)
					break;
			for (Atom a:current.getAdjFrame(4,false)) {
				if(came_from.get(current) != a) {
					came_from.put(a,current);
					stack.push(a);
				}
			}
			steps++;
		}
		
		return (steps==6);
	}
	
	@Override
	public boolean equals(Object obj) { //fast equals, does NOT check changes in bonds etc ...
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof Atom))
			return false;
		Atom other = (Atom) obj;
		if (no != other.no)
			return false;
		if (type != other.type)
			return false;
		return true;
	}

}
