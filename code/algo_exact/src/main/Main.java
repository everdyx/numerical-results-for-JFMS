package main;


public class Main {

	public static void main(String[] args) throws Exception {
		int start_ins; 
		int end_ins;
		int numProd;
		int numMach; 
		int numVeh;
		double cycTime; 
		boolean charge; 
		boolean given_q;
		double timeLimit;
		Main_Exact main_Exact=new Main_Exact();
		
		start_ins=0;
		end_ins=100;
		numProd=2;
		numMach=2;
		numVeh=2;
		cycTime=3;
		charge=true;
		given_q=true;
		timeLimit=15*60;
		main_Exact.run(start_ins, end_ins, numProd, numMach, numVeh, cycTime, charge, given_q,timeLimit);
		System.out.println();
		System.out.println();
		
		cycTime=8;
		main_Exact.run(start_ins, end_ins, numProd, numMach, numVeh, cycTime, charge, given_q,timeLimit);
		System.out.println();
		System.out.println();
		
		cycTime=13;
		main_Exact.run(start_ins, end_ins, numProd, numMach, numVeh, cycTime, charge, given_q,timeLimit);
		System.out.println();
		System.out.println();
		
		numProd=3;
		cycTime=3;
		main_Exact.run(start_ins, end_ins, numProd, numMach, numVeh, cycTime, charge, given_q,timeLimit);
		System.out.println();
		System.out.println();
		
		cycTime=8;
		main_Exact.run(start_ins, end_ins, numProd, numMach, numVeh, cycTime, charge, given_q,timeLimit);
		System.out.println();
		System.out.println();
		
		cycTime=13;
		main_Exact.run(start_ins, end_ins, numProd, numMach, numVeh, cycTime, charge, given_q,timeLimit);
		System.out.println();
		System.out.println();
		
		main_Exact.run(0, 100, 20, 2, 2, 3, charge, given_q,timeLimit);
				
		main_Exact.run(0, 100, 20, 2, 2, 8, charge, given_q,timeLimit);
		
		main_Exact.run(0, 100, 20, 2, 2, 13, charge, given_q,timeLimit);
		
		main_Exact.run(0, 100, 20, 8, 2, 3, charge, given_q,timeLimit);
		
		main_Exact.run(0, 100, 20, 8, 8, 3, charge, given_q,timeLimit);
		
		main_Exact.run(0, 100, 20, 8, 14, 3, charge, given_q,timeLimit);
		
		main_Exact.run(0, 100, 20, 8, 20, 3, charge, given_q,timeLimit);

	}

}
