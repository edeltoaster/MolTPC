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
import java.io.*;

import base.BabelInterface;
import base.Molecule;
import base.RingSearch;

public class IsomerSearch {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Runtime rt = Runtime.getRuntime();
		
		if (args.length==0  || args[0].equals("-help") || args[0].equals("-h")){
			System.out.println("isomersearch [INFOLDER or FILE] STEPS OUTFOLDER");
			System.exit(0);
		}
		
		if (args.length<3){
			System.out.println("wrong arg length");
			System.out.println("isomersearch [INFOLDER or FILE] STEPS OUTFOLDER");
			System.exit(1);
		}
		
		File file = new File(args[0]);
		List<Molecule> molecules = null;
		
		if (file.isDirectory())
			molecules = Molecule.getMoleculesFromFolder(args[0]);
		else {
			molecules = new LinkedList<Molecule>();
			molecules.add(new Molecule(args[0]));
		}
		
		int steps = Integer.parseInt(args[1]);
		
		/*
		 * search
		 */
		
		System.out.println("Isomer-search started for "+args[0]+" with "+steps+" steps");
		System.out.println("Task running on "+BabelInterface.getHostname());
		
		long start=System.currentTimeMillis();
		RingSearch rs = new RingSearch(molecules, steps, 1.0);
		
		System.out.println("Time for searching: "+(System.currentTimeMillis()-start)/1000/60+"min");
		System.out.println(rs.getNoFoundRingClosures()+" ring-conformations found");
		System.out.println(rs.getNoIsomers()+" total isomers");
		System.out.println("Memory used: "+(rt.totalMemory()-rt.freeMemory())/1024/1024+"mb of "+rt.totalMemory()/1024/1024+"mb");
		
		List<Molecule> best_ones = rs.getIsomers();
		rs = null;
		rt.gc();
		
		/*
		 * minimization
		 */
		
		System.out.println("Postprocessing started.");
		best_ones = BabelInterface.batch_cg_minimize(best_ones, 15000, 1);
		Molecule.writeMoleculesToFolder(best_ones,args[2]);
		
		System.out.println("Memory used: "+(rt.totalMemory()-rt.freeMemory())/1024/1024+"mb out of "+rt.totalMemory()/1024/1024);
		System.out.println("Time all: "+(System.currentTimeMillis()-start)/1000/60+"min");
	}

}
