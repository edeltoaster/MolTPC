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

package tools;
import java.util.LinkedList;
import java.util.List;

import base.Molecule;
import base.ProtIt;
import base.TautIt;

public class Enumerate {

	/**
	 * builds isomers -> tautomers + protonation variants
	 * enumerate MOLECULE TAUTOMERISM-FLAG CHARGE OUTPUT-FOLDER
	 * @param args
	 */
	public static void main(String[] args) {
		/**
		 * builds isomers -> tautomers + protonation variants
		 * enumerate MOLECULE TAUTOMERISM-FLAG CHARGE OUTPUT-FOLDER
		 * **/
		if (args.length==0  || args[0].equals("-help") || args[0].equals("-h")){
			System.out.println("buildisomers MOLECULE TAUTOMERISM-FLAG STEREO-FLAG CHARGE OUTFOLDER");
			System.out.println("where taut-flag means:");
			System.out.println("0 : only non-aromatic");
			System.out.println("1 : aromatic with sp2-carbons preserved");
			System.out.println("2 : aromatic without sp2-carbons preserved");
			System.out.println("where stereo-flag means:");
			System.out.println("1 : stereoforms are generated");
			System.out.println("other : no stereoforms generated");
			System.out.println("non-zero CHARGE will invoke the ProtIt-class (+- integers)");
			System.exit(0);
		}
		
		if (args.length < 5){
			System.out.println("wrong arg length!");
			System.out.println("buildisomers MOLECULE TAUTOMERISM-FLAG STEREO-FLAG CHARGE OUTPUT-FOLDER");
			System.out.println("where taut-flag means:");
			System.out.println("0 : only non-aromatic");
			System.out.println("1 : aromatic with sp2-carbons preserved");
			System.out.println("2 : aromatic without sp2-carbons preserved");
			System.out.println("where stereo-flag means:");
			System.out.println("1 : stereoforms are generated");
			System.out.println("other : no stereoforms generated");
			System.out.println("non-zero CHARGE will invoke the ProtIt-class (+- integers)");
			System.exit(1);
		}
		
		Molecule m = new Molecule(args[0]);
		int mode = Integer.parseInt(args[1]);
		int charge = Integer.parseInt(args[3]);
		boolean stereo = args[2].equals("1");
		
		TautIt it = null;
		switch (mode) {
		case 1:
			it = new TautIt(m, true, true, stereo, true, false);
			break;
		case 2:
			it = new TautIt(m, true, false, stereo, true, false);
			break;
		default:
			it = new TautIt(m, false, true, stereo, true, false);
			break;
		}
		
		List<Molecule> isomers = new LinkedList<Molecule>();
		if (charge != 0)
			for (Molecule mol:it.getAllForms()) {
				isomers.addAll(new ProtIt(mol, charge).getAllForms());
			}
		else 
			isomers.addAll(it.getAllForms());
		
		Molecule.writeMoleculesToFolder(isomers, args[4]);
		System.out.println(isomers.size()+" static isomers found.");
		
	}

}
