package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

import basics.Instance;
import ilog.cplex.IloCplex;
import subgradient_method.LogiTaskSortByTss;
import subgradient_method.ProdRatio;
import subgradient_method.ProdSortByWeightProcessTimeRatio;
import subgradient_method.SubgradMethod;
import subgradient_method.TssWithoutLogiCons;

public class Main {

	public static void main(String[] args) throws Exception {
		int start_ins,end_ins,numProd,numVeh;
		
		Main_Large l1=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=10;
		numVeh=2;
		l1.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();
		
		Main_Large l2=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=10;
		numVeh=8;
		l2.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();
		
		Main_Large l3=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=10;
		numVeh=14;
		l3.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();
		
		Main_Large l4=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=10;
		numVeh=20;
		l4.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();

		Main_Large l5=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=20;
		numVeh=2;
		l5.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();
		
		Main_Large l6=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=20;
		numVeh=8;
		l6.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();
		
		Main_Large l7=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=20;
		numVeh=14;
		l7.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();
		
		Main_Large l8=new Main_Large();
		start_ins=0;
		end_ins=100;
		numProd=20;
		numVeh=20;
		l8.run(start_ins,end_ins,numProd,numVeh);
		System.out.println();
		System.out.println();

	}

}
