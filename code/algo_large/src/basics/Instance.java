package basics;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class Instance {
	public int numVertex; // 所有点集合n（包括配送中心和客户点，0为配送中心）
	int numTask;
	public int numVeh; // 车辆数
	public int numProd;
	public int numMach;
	public int bigNum=100000;
	double depot_xcoord=0;
	double depot_ycoord=0;
	public double vehSpeed=1;//m/s
	public double lineSpeed=1;
	public double batteryCapacity=1.6;//1.6;//kWh
	public double carryCapacity=1000;//1000;//kg
	public double iniBattery=0;
	public int numCharVisit;
	public double chargingRate;//=0*7200/1.6; //(s/kWh) //充电量乘以chargingRate=充电时长
	public double consumingRate;//=0*1.6/23040;//(kWh/m)//距离乘以consumingRate=耗电量
//	public boolean considerCharging=false;
	public double jobshopLength=200;
	double jobshopWidth=60;
	public double moldChangeTime=0;//换模时间单位秒
	public Depot depot;
	public ArrayList<Machine> machSet=new ArrayList<Machine>();
	public ArrayList<Product> prodSet=new ArrayList<Product>();
	public ArrayList<LogiTask> logiTaskSet=new ArrayList<LogiTask>();
	public ArrayList<Order> orderSet=new ArrayList<Order>();
	public int[] demands;//=new int[numProd+1];//= {29,28,0}; // 需求量
	public int[][] quanProdMach;//=new int[numProd+1][numMach];//={{30,50},{17,19},
//								 {0,0}};
	public double[] cycTime;//=new double[numProd+1];//= {8,8,0};
	public int[][] K;//=new int[numProd+1][numMach];
	public HashMap<String,LogiTask> PMKtoLogiTask = new HashMap<String,LogiTask>();
	public HashMap<Integer,LogiTask> ItoLogiTask = new HashMap<Integer,LogiTask>();
	public HashMap<Integer,ArrayList<int[]>> PtoArrayofPMKI = new HashMap<Integer,ArrayList<int[]>>();
//	public double moldChangeTime;
	public int maxK=Integer.MIN_VALUE;
	public double objValFactor=1;
	public double cycleTime;
	
	public Instance() {

	}
	
	public Instance(int numProd, int numMach, int numVeh,boolean charge,double vehSpeed,double lineSpeed,double cycleTime, String dataFile) {
		this.numProd=numProd;
		this.numMach=numMach;
		this.numVeh=numVeh;
		this.vehSpeed=vehSpeed;
		this.lineSpeed=lineSpeed;
		this.cycleTime=cycleTime;
		demands=new int[numProd+1];//= {29,28,0}; // 需求量
		quanProdMach=new int[numProd+1][numMach];//={{30,50},{17,19},
//		 {0,0}};
		cycTime=new double[numProd+1];//= {8,8,0};
		K=new int[numProd+1][numMach];
		if (charge) {
			chargingRate=7200/1.6;//((7*3600/12)/vehSpeed)/2;//7200/1.6;
			consumingRate=1.6/23040;//2/(7*3600/vehSpeed);//1.6/23040;
		}
		else {
			chargingRate=0;
			consumingRate=0;
		}
		ReadDataFromFile(dataFile);
	}
	
	private void ReadDataFromFile(String dataFile) {
		try {
//			System.out.println(dataFile);
			Scanner in = new Scanner(new FileReader(dataFile));
			int cnt=0;
			in.nextLine();
			while(in.hasNextLine())
			{	
				Order order=new Order();
				order.id      = cnt;
				order.demand  = in.nextInt();//40~60的均匀分布
//				order.cycleTime    = in.nextInt();
				order.cycleTime=(int) this.cycleTime;
//				order.totalWeight    = in.nextDouble();
				order.weight=new double[numMach];
				for (int i=0;i<numMach;i++) {
					order.weight[i]=in.nextInt();//10+3*(in.nextInt()-15);//demand.txt文件的in.nextInt()是15~25的均匀分布
				}
				in.nextLine();
				orderSet.add(order);
				cnt++;
			}// end for customers
			
			in.close();
			
		} catch (FileNotFoundException e) {
			// File not found
			System.out.println("File not found!");
			System.exit(-1);
		}
	}
	
	public void Initialization(int orderNumStart) {//参数初始化
		//初始化订单
		ArrayList<Integer> pSet =new ArrayList<Integer>();
		for (int i = orderNumStart; i < orderNumStart+numProd; i++) {
//			pSet.add((int)(10000*Math.random()));
			pSet.add(i);
		}
		
		for (int p = 0; p < numProd; p++) {
			demands[p]=orderSet.get(pSet.get(p)).demand;
			cycTime[p]=orderSet.get(pSet.get(p)).cycleTime;
			for (int j = 0; j < numMach; j++) {
				quanProdMach[p][j]=(int)(carryCapacity/orderSet.get(pSet.get(p)).weight[j]);
			}
		}
		//初始化depot
		depot = new Depot();
		depot.setCoord(depot_xcoord, depot_ycoord);
		depot.setId(0);
		//初始化machine
		machSet=new ArrayList<Machine>();
		for (int m = 0; m < numMach; m++) {
			Machine mach = new Machine();
			mach.setCoord((m+1) * jobshopLength / (numMach+1), jobshopWidth / 2);
			mach.setId(m);
			machSet.add(mach);
		}
//		numVertex=numMach+1;
		prodSet=new ArrayList<Product>();
		for (int p = 0; p < numProd; p++) {
			Product prod=new Product();
			prod.setDemand(demands[p]);
			prod.setCycleTime(cycTime[p]);
			prod.setId(p);
			prodSet.add(prod);
		}
		
		logiTaskSet=new ArrayList<LogiTask>();
		int cnt = 0;
		for (Product prod: prodSet) {//遍历P
			for (Machine mach : machSet) {//遍历M
				int numBatch=prod.getDemand() / quanProdMach[prod.getId()][mach.getId()];
				if (numBatch*quanProdMach[prod.getId()][mach.getId()]<prod.getDemand()) {
					numBatch++;
				}
				K[prod.getId()][mach.getId()]=numBatch;
				if (numBatch>maxK) {
					maxK=numBatch;
				}
				for (int k=0;k<numBatch;k++) {//遍历K
					LogiTask logiTask=new LogiTask();
					logiTask.setId(cnt);
					logiTask.setProd(prod);
					logiTask.setMach(mach);
					logiTask.setBatchId(k);
					if (k==numBatch-1) {
						logiTask.setBatchQuant(prod.getDemand()-k*quanProdMach[prod.getId()][mach.getId()]);
					}
					else {
						logiTask.setBatchQuant(quanProdMach[prod.getId()][mach.getId()]);
					}
					logiTask.setAssemTime(logiTask.getBatchQuant()*prod.getCycleTime());
					logiTask.setTransTime((getDistance(depot, logiTask.getMach())+getDistance(logiTask.getMach(),depot))/vehSpeed);
					logiTask.setProcessingTime((getDistance(depot, logiTask.getMach())+getDistance(logiTask.getMach(),depot))/vehSpeed+
							(getDistance(depot, logiTask.getMach())+getDistance(logiTask.getMach(),depot))*consumingRate*chargingRate);
					logiTask.setChargingTime((getDistance(depot, logiTask.getMach())+getDistance(logiTask.getMach(),depot))*consumingRate*chargingRate);
					logiTask.setLambda(0);//暂时设置每个任务的lambda为1,对应的是lagrangian松弛问题中的拉格朗日乘子
					logiTaskSet.add(logiTask);
					cnt++;
				}
			}
		}
		numTask=cnt;
		if (chargingRate>0) {
			numCharVisit=logiTaskSet.size();
		}else {
			numCharVisit=0;
		}
		for (int m = 0; m < numMach; m++) {
			K[numProd][m]=0;
		}
		
		for (LogiTask logiTask:logiTaskSet) {
			
			String PMKstr=logiTask.getProd().getId()+"_"+logiTask.getMach().getId()+"_"+logiTask.getBatchId();
			ItoLogiTask.put(logiTask.getId(),logiTask);
			PMKtoLogiTask.put(PMKstr, logiTask);
			
			PtoArrayofPMKI.put(logiTask.getProd().getId(), new ArrayList<int[]>());
		}
		
		for (LogiTask logiTask:logiTaskSet) {
			int[] PMKI=new int[4];
			PMKI[0]=logiTask.getProd().getId();
			PMKI[1]=logiTask.getMach().getId();
			PMKI[2]=logiTask.getBatchId();
			PMKI[3]=logiTask.getId();
			PtoArrayofPMKI.get(logiTask.getProd().getId()).add(PMKI);
		}
	}
	

	public double getDistance(Vertex a, Vertex b) {
		return Math.abs(a.getXcoord() - b.getXcoord()) + Math.abs(a.getYcoord() - b.getYcoord());
	}


}
