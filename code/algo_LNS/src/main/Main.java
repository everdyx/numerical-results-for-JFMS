package main;


public class Main {

	public static void main(String[] args) throws Exception {

		int numVeh=2;
		int numMach=2;
		int numProd=2;
		double cycTime=3;
		LNS_Large l1=new LNS_Large();
		 l1.run(0,100,true,numProd,numVeh,numMach,cycTime);

		cycTime=8;
		 l1.run(0,100,true,numProd,numVeh,numMach,cycTime);
		cycTime=13;
		l1.run(0,100,true,numProd,numVeh,numMach,cycTime);

		 numProd=3;
		 cycTime=3;
		 l1.run(0,100,true,numProd,numVeh,numMach,cycTime);
		 cycTime=8;
		 l1.run(0,100,true,numProd,numVeh,numMach,cycTime);
		 cycTime=13;
		 l1.run(0,100,true,numProd,numVeh,numMach,cycTime);
		 
		 
		 l1.run(0,100,true,10,2,2,3);
		 
		 l1.run(0,100,true,10,2,2,8);
		 
		 l1.run(0,100,true,10,2,2,13);
		 
		 l1.run(0,100,true,20,2,2,3);
		 
		 l1.run(0,100,true,20,2,2,8);
		 
		 l1.run(0,100,true,20,2,2,13);
		 
		 
		 
		 l1.run(0,100,true,10,2,8,3);
		 
		 l1.run(0,100,true,10,8,8,3);

		 l1.run(0,100,true,10,14,8,3);
		 
		 l1.run(0,100,true,10,20,8,3);
		 
		 l1.run(0,100,true,20,2,8,3);
		 
		 l1.run(0,100,true,20,8,8,3);
		 
		 l1.run(0,100,true,20,14,8,3);
		 
		 l1.run(0,100,true,20,20,8,3);

	}

}
