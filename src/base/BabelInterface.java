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
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.InetAddress;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


public class BabelInterface {
	public static int threads = Runtime.getRuntime().availableProcessors();
	private static String path;
	private static final BabelInterface singleton;
	private static final Runtime rt = Runtime.getRuntime();
	private static String hostname;
	
	static {
		if (path == null || path.equals("")) {
			String[] dirs = {"/usr/local/bin/","/usr/bin/"};
			for (String s:dirs) 
				if (new File(s+"obenergy").exists()) {
					path = s;
					break;
				}
		}
		
		if (path == null) {
			try {
				BufferedReader buffer=new BufferedReader(new FileReader(System.getProperty("user.home")+"/.moltpc/path.cfg"));
				path = buffer.readLine();
				buffer.close();
			} catch (Exception e) {
				System.out.println("No usable OpenBabel-installation found and ~/.moltpc/path.cfg wrong.");
				System.exit(1);
			}
			
			if (!new File(path+"obenergy").exists()) {
				System.out.println("No usable OpenBabel-installation found from line given in ~/.moltpc/path.cfg.");
				System.exit(1);
			}
		}
		
		singleton = new BabelInterface();
		
		try {
			hostname = InetAddress.getLocalHost().getHostName();
		} catch (Exception e) {
			System.out.println("Not able to determine hostname, using 'localhost'.");
			hostname = "localhost";
		}
	}
	
	public static Molecule get_MMFF94_charged(Molecule molecule) {
		Molecule m = null;
		molecule.writeHinFile("smi"+hostname+".hin");
		String[] str = {path+"babel","--partialcharge","mmff94","-i","hin","smi"+hostname+".hin","-o","hin"};
		
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			m = new Molecule(stdout);
			stdout.close();
			pr.waitFor();
			
			// clean
			str = new String[]{"rm","smi"+hostname+".hin"};
			pr = rt.exec(str);
			pr.waitFor();
		} catch (Exception e) {
			System.out.println("Error in Partialcharge-assignment");
			System.out.println(e);
			System.exit(1);
		}
		
		return m;
	}
	
	public static List<Atom> get_chiral_atoms(Molecule molecule) {
		molecule.writeHinFile("smi"+hostname+".hin");
		String[] str = {path+"obchiral","smi"+hostname+".hin"};
		String temp;
		List<Atom> chiral_atoms = new LinkedList<Atom>();
		
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			while ((temp = br.readLine()) != null) {
				if (temp.startsWith("Atom"))
					chiral_atoms.add(molecule.getAtomNo(Integer.parseInt(temp.split("\\s+")[1])));
			}
			
			stdout.close();
			pr.waitFor();
			
		} catch (Exception e) {
			System.out.println("Error in chirality detection");
			System.exit(1);
		}
		
		return chiral_atoms;
	}
	
	public static Molecule sd_minimize(Molecule molecule,int steps,float epsilon) {
		molecule.writeHinFile("mintemp"+hostname+Thread.currentThread().getName()+".hin");
		String[] str = {path+"obminimize","-sd","-ff","mmff94s","-n",Integer.toString(steps),"-o","hin","mintemp"+hostname+Thread.currentThread().getName()+".hin"}; // depends on system
		Molecule minimized = null;
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			
			minimized = new Molecule(stdout);
			stdout.close();
			pr.waitFor();
			
			// clean
			str = new String[]{"rm","mintemp"+hostname+Thread.currentThread().getName()+".hin"};
			pr = rt.exec(str);
			pr.waitFor();
			
		} catch (Exception e) {
			System.out.println("Error in Minimizer-Class");
			System.out.println(e);
			System.exit(1);
		}
		
		return minimized;
	}
	
	public static Molecule cg_minimize(Molecule molecule,int steps,float epsilon) {
		molecule.writeHinFile("mintemp"+hostname+Thread.currentThread().getName()+".hin");
		String[] str = {path+"obminimize","-ff","mmff94s","-n",Integer.toString(steps),"-o","hin","mintemp"+hostname+Thread.currentThread().getName()+".hin"}; // depends on system
		Molecule minimized = null;
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			
			minimized = new Molecule(stdout);
			stdout.close();
			pr.waitFor();
			
			// clean
			str = new String[]{"rm","mintemp"+hostname+Thread.currentThread().getName()+".hin"};
			pr = rt.exec(str);
			pr.waitFor();
			
		} catch (Exception e) {
			System.out.println("Error in Minimizer-Class");
			System.out.println(e);
			System.exit(1);
		}
		
		return minimized;
	}
	
	public static double energy_mmff94(Molecule molecule, float epsilon) {
		molecule.writeHinFile("temp"+hostname+".hin");
		String[] str = { path+"obenergy","-ff","mmff94s", "temp"+hostname+".hin" };
		double energy = 0;
		try {
		
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			String line;
			String temp = null;
			
			// clean up if any output in stdout
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			while ((line = br.readLine()) != null)
				temp = line;
			br.close();
			stdout.close();
			pr.waitFor();
			
			// happens in few cases, why? no clue! concurrency/locking/network! problem? workaround:
			if (temp == null) {
				rt.gc();
				return energy_mmff94(molecule, epsilon);
			}
			energy = Double.parseDouble(temp);
			
			// clean
			str = new String[]{"rm","temp"+hostname+".hin"};
			pr = rt.exec(str);
			pr.waitFor();
			pr.getInputStream().close();
			pr.getOutputStream().close();
			pr.getErrorStream().close();
			
		} catch (Exception e) {
			System.out.println("Error while using obenergy");
			e.printStackTrace();
			System.exit(1);
		}
		
		return energy; //  kcal/mol
	}
	
	public static double energy_mmff94fast(Molecule molecule,float epsilon) {
		molecule.writeHinFile("temp"+hostname+Thread.currentThread().getName()+".hin");
		String[] str = { path+"obenergy","-ff","mmff94s", "temp"+hostname+Thread.currentThread().getName()+".hin" };
		double energy = 0;
		try {
		
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			String line;
			String temp = null;
			
			// clean up if any output in stdout
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			while ((line = br.readLine()) != null) {
				temp = line;
			}
			br.close();
			stdout.close();
			pr.waitFor();
			
			
			// happens in few cases, why? no clue! concurrency/locking/network! problem? workaround:
			if (temp == null) {
				rt.gc();
				return energy_mmff94fast(molecule, epsilon);
			}
			energy = Double.parseDouble(temp);
			
			
		} catch (Exception e) {
			System.out.println("Error while using obenergy");
			e.printStackTrace();
			System.exit(1);
		}
		
		return energy; //  kcal/mol
	}
	
	/**
	 * 
	 * @param molecule
	 * @param epsilon
	 * @return total, bond, angle, torsional, oop, vdw ,estat
	 */
	public static double[] energy_GAFFfast(Molecule molecule,float epsilon) {
		molecule.writeHinFile("temp"+hostname+Thread.currentThread().getName()+".hin");
		String[] str = { path+"obenergy","-ff","gaff","-v", "temp"+hostname+Thread.currentThread().getName()+".hin" };
		double energy = 0;
		double bond_stretch = 0;
		double angle_bend = 0;
		double torsional = 0;
		double inpr_torsional = 0;
		double vdw = 0;
		double estat = 0;
		try {
			
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			String line;
			String[] temp = null;
			
			// clean up if any output in stdout
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			while ((line = br.readLine()) != null) {
				temp = line.split("\\s+");
				if (temp.length>5 && temp[1].equals("TOTAL")){
					if (temp[2].startsWith("BO"))
						bond_stretch = Double.parseDouble(temp[6])*0.2388;
					else if (temp[2].startsWith("AN"))
						angle_bend = Double.parseDouble(temp[6])*0.2388;
					else if (temp[2].startsWith("TO"))
						torsional = Double.parseDouble(temp[5])*0.2388;
					else if (temp[2].startsWith("IM"))
						inpr_torsional = Double.parseDouble(temp[5])*0.2388;
					else if (temp[2].startsWith("VAN"))
						vdw = Double.parseDouble(temp[7])*0.2388;
					else if (temp[2].startsWith("ELEC"))
						estat = Double.parseDouble(temp[5])*0.2388/epsilon;
				}
				
			}
			br.close();
			stdout.close();
			pr.waitFor();
			
			
			// happens in few cases, why? no clue! concurrency/locking/network! problem? workaround:
			if (temp == null) {
				rt.gc();
				return energy_GAFFfast(molecule, epsilon);
			}
			energy = Double.parseDouble(temp[0])*0.2388;// WRONG, take SUM !!!
			
			
		} catch (Exception e) {
			System.out.println("Error while using GAFF-obenergy");
			e.printStackTrace();
			System.exit(1);
		}
		
		energy = bond_stretch+angle_bend+torsional+inpr_torsional+vdw+estat;
		
		return new double[]{energy,bond_stretch,angle_bend,torsional,inpr_torsional,vdw,estat}; //  kcal/mol
	}
	
	/**
	 * 
	 * @param molecule
	 * @param epsilon
	 * @return total, bond, angle, stretch_bend, torsional, oop, vdw ,estat
	 */
	public static double[] energy_MMFF94sCompfast(Molecule molecule,float epsilon) {
		molecule.writeHinFile("temp"+hostname+Thread.currentThread().getName()+".hin");
		String[] str = { path+"obenergy","-ff","mmff94s","-v", "temp"+hostname+Thread.currentThread().getName()+".hin" };
		double energy = 0;
		double bond_stretch = 0;
		double angle_bend = 0;
		double stretch_bend = 0;
		double torsional = 0;
		double oop = 0;
		double vdw = 0;
		double estat = 0;
		try {
			
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			String line;
			String[] temp = null;
			
			// clean up if any output in stdout
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			while ((line = br.readLine()) != null) {
				temp = line.split("\\s+");
				if (temp.length>5 && temp[1].equals("TOTAL")){
					if (temp[2].startsWith("BO"))
						bond_stretch = Double.parseDouble(temp[6]);
					else if (temp[2].startsWith("AN"))
						angle_bend = Double.parseDouble(temp[6]);
					else if (temp[2].startsWith("STR"))
						stretch_bend = Double.parseDouble(temp[6]);
					else if (temp[2].startsWith("TO"))
						torsional = Double.parseDouble(temp[5]);
					else if (temp[2].startsWith("OU"))
						oop = Double.parseDouble(temp[6]);
					else if (temp[2].startsWith("VAN"))
						vdw = Double.parseDouble(temp[7]);
					else if (temp[2].startsWith("ELEC"))
						estat = Double.parseDouble(temp[5]);
				}
				
			}
			br.close();
			stdout.close();
			pr.waitFor();
			
			
			// happens in few cases, why? no clue! concurrency/locking/network! problem? workaround:
			if (temp == null) {
				rt.gc();
				return energy_MMFF94sCompfast(molecule, epsilon);
			}
			energy = Double.parseDouble(temp[0]);
			
			
		} catch (Exception e) {
			System.out.println("Error while using MMFF-obenergy");
			e.printStackTrace();
			System.exit(1);
		}
		
		return new double[]{energy,bond_stretch,angle_bend,stretch_bend,torsional,oop,vdw,estat}; //  kcal/mol
	}
	
	public static String getSMILES(Molecule molecule) {
		molecule.writeHinFile("smi"+hostname+".hin");
		String[] str = {path+"babel","-i","hin","smi"+hostname+".hin","-o","can"}; // depends on system
		String smiles = null;
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			smiles = br.readLine().split("\\s+")[0];
			br.close();
			stdout.close();
			pr.waitFor();
						
			// clean
			str = new String[]{"rm","smi"+hostname+".hin"};
			pr = rt.exec(str);
			pr.waitFor();
		} catch (Exception e) {
			System.out.println("Error while using SMILES");
			e.printStackTrace();
			System.exit(1);
		}
		
		return smiles;
	}
	
	public static String getSMILESfast(Molecule molecule) {
		molecule.writeHinFile("smi"+hostname+Thread.currentThread().getName()+".hin");
		String[] str = {path+"babel","-i","hin","smi"+hostname+Thread.currentThread().getName()+".hin","-o","can"}; // depends on system
		String smiles = null;
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			smiles = br.readLine().split("\\s+")[0];
			br.close();
			stdout.close();
			pr.waitFor();

			// no cleaning!
		} catch (Exception e) {
			System.out.println("Error while using SMILESfast");
			e.printStackTrace();
			System.exit(1);
		}
		
		return smiles;
	}
	
	public static String getINCHI(Molecule molecule) {
		molecule.writeHinFile("smi"+hostname+".hin");
		String[] str = {path+"babel","-i","hin","smi"+hostname+".hin","-o","inchi"}; // depends on system
		String inchi_str = null;
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			inchi_str = br.readLine().split("\\s+")[0];
			br.close();
			stdout.close();
			pr.waitFor();

			// clean
			str = new String[]{"rm","smi"+hostname+".hin"};
			pr = rt.exec(str);
			pr.waitFor();
			
		} catch (Exception e) {
			System.out.println("Error while using INCHI");
			e.printStackTrace();
			System.exit(1);
		}
		
		return inchi_str;
	}
	
	public static String getINCHIfast(Molecule molecule) {
		molecule.writeHinFile("smi"+hostname+Thread.currentThread().getName()+".hin");
		String[] str = {path+"babel","-i","hin","smi"+hostname+Thread.currentThread().getName()+".hin","-o","inchi"}; // depends on system
		String inchi_str = null;
		try {
			Process pr = rt.exec(str);
			InputStream stdout = pr.getInputStream();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
			inchi_str = br.readLine().split("\\s+")[0];
			br.close();
			stdout.close();
			pr.waitFor();
			
			// no cleaning!
		} catch (Exception e) {
			System.out.println("Error while using INCHIfast");
			e.printStackTrace();
			System.exit(1);
		}
		
		return inchi_str;
	}
	
	public static List<Molecule> batch_sd_minimize(Collection<Molecule> molecules,int steps,float epsilon) {
		List<Molecule> minimized = new LinkedList<Molecule>();
		List<SD_Minimizer> tasks = new LinkedList<SD_Minimizer>();
		
		for (Molecule m:molecules)
			tasks.add(singleton.new SD_Minimizer(m, steps, epsilon));
		
		ExecutorService es = Executors.newFixedThreadPool(threads);
		try {
			List<Future<Molecule>> results = es.invokeAll(tasks);
			for (Future<Molecule> r:results)
				minimized.add(r.get());
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		es.shutdown();
		
		minimizeClean();
		
		return minimized;
	}
	
	public static List<Molecule> batch_cg_minimize(Collection<Molecule> molecules,int steps,float epsilon) {
		List<Molecule> minimized = new LinkedList<Molecule>();
		List<CG_Minimizer> tasks = new LinkedList<CG_Minimizer>();
		
		for (Molecule m:molecules)
			tasks.add(singleton.new CG_Minimizer(m, steps, epsilon));
		
		ExecutorService es = Executors.newFixedThreadPool(threads);
		try {
			List<Future<Molecule>> results = es.invokeAll(tasks);
			for (Future<Molecule> r:results)
				minimized.add(r.get());
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		es.shutdown();
		
		minimizeClean();
		
		return minimized;
	}
	
	public static List<String> batch_SMILES(Collection<Molecule> molecules) {
		List<String> minimized = new LinkedList<String>();
		List<SMILES> tasks = new LinkedList<SMILES>();
		
		for (Molecule m:molecules)
			tasks.add(singleton.new SMILES(m));
		
		ExecutorService es = Executors.newFixedThreadPool(threads);
		try {
			List<Future<String>> results = es.invokeAll(tasks);
			for (Future<String> r:results)
				minimized.add(r.get());
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		es.shutdown();
		
		smiClean();
		
		return minimized;
	}
	
	public static Map<Double,Molecule> batch_energy_sorted(Collection<Molecule> molecules,float epsilon) {
		Map<Double,Molecule> sorted_map = new TreeMap<Double, Molecule>();
		List<Energy_Calculator> tasks = new LinkedList<Energy_Calculator>();
		
		for (Molecule m:molecules)
			tasks.add(singleton.new Energy_Calculator(m, epsilon));
		
		ExecutorService es = Executors.newFixedThreadPool(threads);
		try {
			List<Future<Object[]>> results = es.invokeAll(tasks);
			for (Future<Object[]> r:results) {
				Object[] pair = r.get();
				double temp = (Double) pair[1];
				Molecule m = (Molecule) pair[0];
				while (sorted_map.containsKey(temp)) // ensure nothing gets lost
					temp += 0.00000000000001;
					
				sorted_map.put(temp,m);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		es.shutdown();
		
		energyClean();
		
		return sorted_map;
	}
	
	/**
	 *  Batch energy calculation, since ordering of sets not known, pairs of molecule/energy(double)
	 *  are returned in an Object-tuple
	 * @param molecules
	 * @param epsilon
	 * @return
	 */
	public static List<Object[]> batch_energy(Collection<Molecule> molecules,float epsilon) {
		List<Energy_Calculator> tasks = new LinkedList<Energy_Calculator>();
		List<Object[]> mol_energies = new LinkedList<Object[]>();
		
		for (Molecule m:molecules)
			tasks.add(singleton.new Energy_Calculator(m, epsilon));
		
		ExecutorService es = Executors.newFixedThreadPool(threads);
		
		try {
			List<Future<Object[]>> results = es.invokeAll(tasks);
			
			for (Future<Object[]> r:results) 
				mol_energies.add(r.get());
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		es.shutdown();
		
		energyClean();
		
		return mol_energies;
	}
	
	class SD_Minimizer implements Callable<Molecule> {
		private Molecule molecule;
		private int steps;
		private float epsilon;
		
		public SD_Minimizer(Molecule molecule,int steps,float epsilon) {
			this.molecule = molecule;
			this.steps = steps;
			this.epsilon = epsilon;
		}
		
		@Override
		public Molecule call() {
			molecule.writeHinFile("mintemp"+hostname+Thread.currentThread().getName()+".hin");
			String name = molecule.getName();
			String[] str = {path+"obminimize","-sd","-ff","mmff94s","-n",Integer.toString(steps),"-o","hin","mintemp"+hostname+Thread.currentThread().getName()+".hin"}; // depends on system
			
			try {
				Process pr = rt.exec(str);
				InputStream stdout = pr.getInputStream();
				molecule = new Molecule(stdout);
				stdout.close();
				pr.waitFor();
				
			} catch (Exception e) {
				System.out.println("Error in Minimizer-Class");
				System.out.println(e);
				System.exit(1);
			}
			
			molecule.setName(name);
			return molecule;
		}
		
	}
	
	class CG_Minimizer implements Callable<Molecule> {
		private Molecule molecule;
		private int steps;
		private float epsilon;
		
		public CG_Minimizer(Molecule molecule,int steps,float epsilon) {
			this.molecule = molecule;
			this.steps = steps;
			this.epsilon = epsilon;
		}
		
		@Override
		public Molecule call() {
			molecule.writeHinFile("mintemp"+hostname+Thread.currentThread().getName()+".hin");
			String name = molecule.getName();
			String[] str = {path+"obminimize","-ff","mmff94s","-n",Integer.toString(steps),"-o","hin","mintemp"+hostname+Thread.currentThread().getName()+".hin"}; // depends on system
			
			try {
				Process pr = rt.exec(str);
				InputStream stdout = pr.getInputStream();
				molecule = new Molecule(stdout);
				stdout.close();
				pr.waitFor();
				
			} catch (Exception e) {
				System.out.println("Error in Minimizer-Class");
				System.out.println(e);
				System.exit(1);
			}
			
			molecule.setName(name);
			return molecule;	
		}
		
	}
	
	class Energy_Calculator implements Callable<Object[]> {
		
		private Molecule molecule;
		private float epsilon;
		
		public Energy_Calculator(Molecule molecule, float epsilon) {
			this.molecule = molecule;
			this.epsilon = epsilon;
		}
		
		@Override
		public Object[] call() throws Exception {
			String temp = null;
			double energy = 0;
			molecule.writeHinFile("temp"+hostname+Thread.currentThread().getName()+".hin");
			String[] str = { path+"obenergy","-ff","mmff94s", "temp"+hostname+Thread.currentThread().getName()+".hin" };
			
			// happens in few cases, why? no clue! concurrency/locking/network! problem? workaround:
			while (temp == null) {
				
				try {
					Process pr = rt.exec(str);
					InputStream stdout = pr.getInputStream();
					
					
					// clean up if any output in stdout
					BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
					
					temp = br.readLine();
					
					br.close();
					stdout.close();
					pr.waitFor();
					
					energy = Double.parseDouble(temp);
					
				} catch (Exception e) {
					temp = null; // if it was another exception ..
					Thread.sleep(500);
				}
			}
			
			
			// build "single entry map"
			Object[] pair = new Object[2];
			pair[0] = molecule;
			pair[1] = energy;
			
			return pair; //  kcal/mol
		}
		
	}
	
	class SMILES implements Callable<String> {
		private Molecule molecule;
		
		
		public SMILES(Molecule molecule) {
			this.molecule = molecule;
		}
		
		@Override
		public String call() {
			molecule.writeHinFile("smi"+hostname+Thread.currentThread().getName()+".hin");
			String[] str = {path+"babel","-i","hin","smi"+hostname+Thread.currentThread().getName()+".hin","-o","can"}; // depends on system
			String smiles = null;
			try {
				Process pr = rt.exec(str);
				InputStream stdout = pr.getInputStream();
				
				BufferedReader br = new BufferedReader(new InputStreamReader(stdout));
				smiles = br.readLine().split("\\s+")[0];
				br.close();
				stdout.close();
				pr.waitFor();

				// no cleaning!
			} catch (Exception e) {
				System.out.println("Error while using batchSMILESfast");
				e.printStackTrace();
				System.exit(1);
			}
			
			return smiles;	
		}
		
	}
	
	public static void minimizeClean() {
		String[] str = {"/bin/bash", "-c","rm -rf mintemp"+hostname+"*"};
		
		// clean
		Process pr;
		try {
			pr = rt.exec(str);
			pr.waitFor();
			pr.getInputStream().close();
			pr.getOutputStream().close();
			pr.getErrorStream().close();
		} catch (Exception e) {
			System.out.println("Error in Minimizer-Cleaner");
			e.printStackTrace();
			System.exit(1);
		}
		
	}

	public static void energyClean() {
		String[] str = {"/bin/bash", "-c","rm -rf temp"+hostname+"*"};
		
		// clean
		Process pr;
		try {
			pr = rt.exec(str);
			pr.waitFor();
			pr.getInputStream().close();
			pr.getOutputStream().close();
			pr.getErrorStream().close();
		} catch (Exception e) {
			System.out.println("Error in Energy-Cleaner");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static void smiClean() {
		String[] str = {"/bin/bash", "-c","rm -rf smi"+hostname+"*"};
		
		// clean
		Process pr;
		try {
			pr = rt.exec(str);
			pr.waitFor();
			pr.getInputStream().close();
			pr.getOutputStream().close();
			pr.getErrorStream().close();
		} catch (Exception e) {
			System.out.println("Error in Smi-Cleaner");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static String getHostname() {
		return hostname;
	}
	
}
