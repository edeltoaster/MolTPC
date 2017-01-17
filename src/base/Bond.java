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
import java.util.List;


public class Bond {
	private Atom a1,a2;
	private byte type;
	private boolean marked;
	
	public Bond(Atom a1,Atom a2,byte type) {
		this.a1 = a1;
		this.a2 = a2;
		this.type = type;
		a1.addBond(this);
		a2.addBond(this);
	}
	
	/**
	 * Returns the atom connected by this bond to the given atom
	 * 
	 * @param atom atom known
	 * @return the connected atom
	 */
	public Atom getOtherAtom(Atom atom) {
		if (atom.equals(a1))
			return a2;
		else 
			return a1;
	}
	
	public byte getType() {
		return type;
	}
	
	public void setType(byte type) {
		this.type = type;
	}
	
	public void setType(int type) {
		this.type = (byte) type;
	}
	
	public void mark() {
		this.marked = true;
	}
	
	public void unmark() {
		this.marked = false;
	}
	
	public boolean isMarked() {
		return this.marked;
	}
	
	public void delBond() {
		a1.deleteBond(this);
		a2.deleteBond(this);
	}
	
	public void setBindToInsteadOf(Atom to,Atom not_anymore) {
		not_anymore.deleteBond(this);
		if (not_anymore.equals(a1)) {
			a1 = to;
			to.addBond(this);
			}
		else {
			a2 = to;
			to.addBond(this);
		}
	}
	
	public String toString() {
		return ("B:"+a1.getNo()+" to "+a2.getNo()+" with "+type);
	}
	
	public String getStrType() {
		switch (type) {
		case 1:
			return "s";
		case 2:
			return "d";
		case 3:
			return "t";
		default:
			return "a";
		}
	}
	
	public Atom getA1() {
		return a1;
	}

	public Atom getA2() {
		return a2;
	}

	public void flip() {
		
		switch (type) {
		case 4: // aromatic bonds
			break;
		case 1:
			type = 2;
			break;
		default: // db/tb
			type = (byte) (type - 1);
			break;
		}
	}
	
	public List<Atom> getMinimalRotAtoms() {
		List<Atom> a1_reach = a1.getRotAtoms(this);
		List<Atom> a2_reach = a2.getRotAtoms(this);
		
		if (a1_reach.size()<a2_reach.size())
			return a1_reach;
		
		return a2_reach;
	}
	
	public double getAngle() {
		float[] vec1 = new float[3];
		float[] vec2 = new float[3];
		float[] vec3 = new float[3];
		float x1,x2,x3,y1,y2,y3;
		
		// set neighbour of a1 and a2
		Atom a = null;
		Atom b = null;
		
		for (Atom temp:a1.getAdjAtoms()){
			if (!temp.equals(a2) && temp.getType()==Element.C) {
				a = temp;
				break;
			}
		}
		if (a==null)
			a = a1.getAdjAtoms().get(0);
		if (a.equals(a2))
			a = a1.getAdjAtoms().get(1);
		
		
		for (Atom temp:a2.getAdjAtoms()){
			if (!temp.equals(a1) && temp.getType()==Element.C) {
				b = temp;
				break;
			}
		}
		if (b==null)
			b = a2.getAdjAtoms().get(0);
		if (a.equals(a1))
			a = a2.getAdjAtoms().get(1);
		
		// normalvector of first plane
		x1 = a.getX() - a1.getX(); x2 = a.getY() - a1.getY(); x3 = a.getZ() - a1.getZ();
		y1 = a1.getX() - a2.getX(); y2 = a1.getY() - a2.getY(); y3 = a1.getZ() - a2.getZ();
		vec1[0] = x2*y3 - x3*y2;
		vec1[1] = x3*y1 - x1*y3;
		vec1[2] = x1*y2 - x2*y1;

		// normalvector of second plane
		x1 = a1.getX() - a2.getX(); x2 = a1.getY() - a2.getY(); x3 = a1.getZ() - a2.getZ();
		y1 = a2.getX() - b.getX(); y2 = a2.getY() - b.getY(); y3 = a2.getZ() - b.getZ();
		vec2[0] = x2*y3 - x3*y2;
		vec2[1] = x3*y1 - x1*y3;
		vec2[2] = x1*y2 - x2*y1;

		// calculate angle between these vectors
		double angle = Math.toDegrees(Math.acos((vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]) / 
				(Math.sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])
						*Math.sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]))));
			
		vec3[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
		vec3[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
		vec3[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
		
		return Math.signum(((a1.getX() - a2.getX())*vec3[0]
				+a1.getY() - a2.getY()*vec3[1]+
				a1.getZ() - a2.getZ()*vec3[2]))*(-1)*angle;
	}
	
	public int[] getAngleAtoms() {
		
		// set neighbour of a1 and a2
				Atom a = null;
				Atom b = null;
				
				for (Atom temp:a1.getAdjAtoms()){
					if (!temp.equals(a2) && temp.getType()==Element.C) {
						a = temp;
						break;
					}
				}
				if (a==null)
					a = a1.getAdjAtoms().get(0);
				if (a.equals(a2))
					a = a1.getAdjAtoms().get(1);
				
				
				for (Atom temp:a2.getAdjAtoms()){
					if (!temp.equals(a1) && temp.getType()==Element.C) {
						b = temp;
						break;
					}
				}
				if (b==null)
					b = a2.getAdjAtoms().get(0);
				if (a.equals(a1))
					a = a2.getAdjAtoms().get(1);
				
		int[] atoms = new int[]{a.getNo(),a1.getNo(),a2.getNo(),b.getNo()};
		
		return atoms;
	}
	
	public void rotate(float degrees) {
		float x = a2.getX()-a1.getX();
		float y = a2.getY()-a1.getY();
		float z = a2.getZ()-a1.getZ();
		
		// for unit-quaternion
		float factor = (float) (1/Math.sqrt(x*x+y*y+z*z));
		
		float[] rot_vector = new float[]{factor*x, factor*y, factor*z};
		
		// the quaternion rotates with "step"-width in every iteration
		float sin_coeff = (float) Math.sin((degrees*Math.PI)/360);
		float cos_coeff = (float) Math.cos((degrees*Math.PI)/360);
		
		Quaternion rot_quat = new Quaternion(cos_coeff,
				sin_coeff*rot_vector[0], 
				sin_coeff*rot_vector[1],
				sin_coeff*rot_vector[2]);
		
		Quaternion rot_quat_conj = rot_quat.getConjQuat();
		
		// the "correction"-vector
		float[] corr = new float[]{a1.getX(),a1.getY(),a1.getZ()};
		
		for (Atom a:getMinimalRotAtoms())
			a.rotate(rot_quat, rot_quat_conj, corr);
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((a1 == null) ? 0 : a1.getNo());
		result = prime * result + ((a2 == null) ? 0 : a2.getNo());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof Bond))
			return false;
		Bond other = (Bond) obj;
		if (type!=other.type)
			return false;
		if (a1 == null) {
			if (other.a1 != null)
				return false;
		} else if (!a1.equals(other.a1))
			return false;
		if (a2 == null) {
			if (other.a2 != null)
				return false;
		} else if (!a2.equals(other.a2))
			return false;
		
		return true;
	}

}
