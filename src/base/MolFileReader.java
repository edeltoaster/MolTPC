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
import java.util.LinkedList;
import java.util.List;

public abstract class MolFileReader {
	
	protected Atom[] frame;
	protected List<Atom> protons=new LinkedList<Atom>();
	protected String name;
	protected boolean explHs;
	
	public Atom[] getFrame() {
		return frame;
	}

	public List<Atom> getProtons() {
		return protons;
	}

	public String getName() {
		return name;
	}

	public boolean isExplHs() {
		return explHs;
	}
	
}
