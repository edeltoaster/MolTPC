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

public class Quaternion {
	private float a;
	private float b;
	private float c;
	private float d;
	
	public Quaternion(float a,float b,float c,float d) {
		this.a = a;
		this.b = b;
		this.c = c;
		this.d = d;
	}
	
	// for points -> Quaternion
	public Quaternion(float x,float y,float z) {
		this.a = 0;
		this.b = x;
		this.c = y;
		this.d = z;
	}
	
	public Quaternion getConjQuat() {
		return new Quaternion(a,-b,-c,-d);
	}
	
	public Quaternion mult(Quaternion q2) {
		return new Quaternion((a*q2.a - b*q2.b - c*q2.c - d*q2.d), 
				(a*q2.b + b*q2.a + c*q2.d - d*q2.c), 
				(a*q2.c - b*q2.d + c*q2.a + d*q2.b), 
				(a*q2.d + b*q2.c - c*q2.b + d*q2.a));
	}
	
	public float[] point_mult(Quaternion q2) {
		float[] point = new float[3];
		// x
		point[0] = a*q2.b + b*q2.a + c*q2.d - d*q2.c;
		// y
		point[1] = a*q2.c - b*q2.d + c*q2.a + d*q2.b;
		// z
		point[2] = a*q2.d + b*q2.c - c*q2.b + d*q2.a;
		
		return point;
	}
	public float getA() {
		return a;
	}

	public void setA(float a) {
		this.a = a;
	}

	public float getB() {
		return b;
	}

	public void setB(float b) {
		this.b = b;
	}

	public float getC() {
		return c;
	}

	public void setC(float c) {
		this.c = c;
	}

	public float getD() {
		return d;
	}

	public void setD(float d) {
		this.d = d;
	}
	
}
