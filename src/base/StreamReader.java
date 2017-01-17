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
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;


public class StreamReader extends MolFileReader {

	/**
	 * Reads a .hin file and delivers getters for Molecule-attributes
	 * @param file path of the .hin-file
	 */
	public StreamReader(InputStream input) {
		List<Atom> all_atoms = new LinkedList<Atom>();
		List<String> bond_info = new LinkedList<String>();
		
		try {
			BufferedReader buffer=new BufferedReader(new InputStreamReader(input));
			String line;
			String[] temp;
			
			while ((line = buffer.readLine()) != null) {
				
				if (line.equals(";")) //comments
					continue;
				if (line.startsWith("endmol")) //files with many molecules?
					break;
				
				temp=line.split("\\s+");
				Atom temp_atom;
				if (temp[0].equals("atom")) {
					temp_atom = new Atom(Element.valueOf(temp[3]),Integer.parseInt(temp[1]),Float.parseFloat(temp[7]),Float.parseFloat(temp[8]),Float.parseFloat(temp[9]),Float.parseFloat(temp[6]));
					String bonds = "";
					for (int i = 11;(i-11) < Integer.parseInt(temp[10])*2 ;i += 2){
						
						// matches all_atoms index with -1
						bonds += (Integer.parseInt(temp[i])-1)+" ";
						
						
						if (temp[i+1].equals("s"))
							bonds += "1";
						else if (temp[i+1].equals("d"))
							bonds += "2";
						else if (temp[i+1].equals("t"))
							bonds += "3";
						else
							bonds += "4"; //aromatic
						
						//separate bonds
						bonds += "\n";
					}
					
					bond_info.add(bonds);
					all_atoms.add(temp_atom);
					
				}
				
				// more?
				// multi-molecules in one hin etc
				
			}
			
			buffer.close();
			
		
		} catch (Exception e) {
			System.out.println("Stream not in valid format.");
			System.exit(1);
		}
		

		//build molecule object stuff out of atoms and bondstring
		for (int i=0;i<bond_info.size();i++) {
			String[] bonds=bond_info.get(i).split("\n");
			for (String s:bonds) {
				int atom2=Integer.parseInt(s.split(" ")[0]);
				byte type=Byte.parseByte(s.split(" ")[1]);
				
				//every bond only once, constructor adds them to both atoms
				if (i<atom2)
					new Bond(all_atoms.get(i),all_atoms.get(atom2),type);
			}
		}
		
		//separate framework/protons
		List<Atom> tempframe=new LinkedList<Atom>();
		for (Atom a:all_atoms)
			if (a.getType()==Element.H)
				this.protons.add(a);
			else
				tempframe.add(a);
		
		//set frame of framework-atoms etc
		this.frame = tempframe.toArray(new Atom[tempframe.size()]);
		this.explHs = (protons.size()>0);
		
		//filename sets name
		this.name = "";
	}
	
}
