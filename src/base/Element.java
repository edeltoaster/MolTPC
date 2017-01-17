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


//Values from GAFF / AMBER
public enum Element {
	H(1,Integer.MAX_VALUE,0.6f,0.0157f),  // DEFAULT Hx
	C(4,0,1.908f,0.086f), // DEFAULT: all but sp3 C
	N(3,1,1.824f,0.17f), //NH2-X will get a score of 2 for protonation
	O(2,-1,1.6612f,0.21f), //COOH OH -> O will get protonation score of -2
	F(1,0,1.75f,0.061f),
	P(3,0,2.1f,0.2f),
	S(2,0,2.0f,0.25f),
	Cl(1,0,1.948f,0.265f),
	Br(1,0,2.22f,0.32f),
	I(1,0,2.35f,0.4f);

	private final int valency;
	private final int loseHEnergy;
	private final float vdW;
	private final float eps;
	
	Element(int val,int loseH,float vdW,float eps) {
		this.valency = val;
		this.loseHEnergy = loseH;
		this.vdW = vdW;
		this.eps = eps;
	}
	
	public int getValency() {
		return valency;
	}
	
	public int getLoseHEnergy() {
		return loseHEnergy;
	}
	
	public float getVdW(Atom a) {
		switch (this) {
			case H:
				Atom adj = a.getAdjAtoms().get(0);
				if (adj.getType()==Element.C) { // many types of Hs
					if (adj.isAromatic())
						return 1.459f; // ha
					
					int elwd = 0; 
					for (Atom b:adj.getAdjAtoms()) // counting electrowithdrawing groups
						if (b.isElWithdraw())
							elwd++;
					
					if (adj.getAdjAtoms().size()==4) { // sp3 C
						switch (elwd) {
						case 1:	// h4
							return 1.409f;
						default: // h5
							return 1.359f;
						}
					} else { // otherwise
						switch (elwd) {
						case 0:	 // hc
							return 1.487f;
						case 1:	// h1
							return 1.387f;
						case 2:	// h2
							return 1.287f;
						default: // h3
							return 1.187f;
						}
					}
				}
				else if (adj.getType()==Element.O)
					return 0; // normally zero, eps is 0 for energy-scoring ... 
				
				break;
			case O: // default O=
				if (a.getAdjAtoms().size()==2) {
					if (a.getAdjHs().size()!=0) 
						return 1.721f; // OH
					else 
						return 1.6837f; // OS,ether
				}
				break;
			
			default:
				break;
		}
		return vdW;
	}
	
	public float getEps(Atom a) {
		switch (this) {
			case H:
				Atom adj = a.getAdjAtoms().get(0);
				if (adj.getType()==Element.C) { // many types of Hs
					if (adj.getAdjAtoms().size()<4 && !adj.isAromatic()) // sp3 -> default
						break;
					
					return 0.015f;
					
				}
				else if (adj.getType()==Element.O)
					return 0; // HO (probably manually?)
				
				break;
			case O:
				if (a.getAdjAtoms().size()==2) {
					if (a.getAdjHs().size()!=0) 
						return 0.2104f; // OH
					else 
						return 0.17f; // OS,ether
				}
				break;
			case C: // default: all but sp3
				if (a.getAdjAtoms().size()==4)
					return 0.1094f;
				break;
			default:
				break;
		}
		return eps;
	}
}
