package main;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.io.FileWriter;
import basics.Instance;
import subgradient_method.LogiTaskSortByTss;
import subgradient_method.ProdRatio;
import subgradient_method.ProdSortByWeightProcessTimeRatio;
import subgradient_method.SubgradMethod;
import subgradient_method.TssWithoutLogiCons;

public class LNS_Large {
	public void run(int start_ins, int end_ins, boolean onlySeq,int numProd,int numVeh,int numMach,double cycTime) throws ClassNotFoundException, IOException{
		// boolean onlySeq = true;
		// int start_ins = 1;
		// int end_ins = start_ins + 100;
		boolean charge;
		// int numProd = 10;
		// double cycTime = 3;
		double vehSpeed = 1;// m/s

		if (numProd <= 3) {
			charge = true;
		} else {
			charge = true;
		}

		// int numMach = 8;
		// int numVeh = 2;
		double lineSpeed = 0.5;
		String dataFile;
		if (System.getProperty("os.name").toLowerCase().contains("win")) {
			dataFile="data\\demand3.txt";
		}
		else {
			dataFile="data/demand3.txt";
		}
		String csvFilePath = ".\\output\\"+"P" + numProd + "M" + numMach + "V" + numVeh+" CT"+cycTime+"_onlySeq_"+onlySeq+".csv";
		File file = new File(csvFilePath);
		if (!file.exists()) {
			boolean fileCreated = file.createNewFile();
			if (fileCreated) {
				System.out.println("File created successfully.");
			} else {
				System.out.println("File already exists or cannot be created.");
			}
		} else {
			System.out.println("File already exists.");
		}
		CSVProcess csvProcess=new CSVProcess();
		System.out.println("P" + numProd + "M" + numMach + "V" + numVeh+"onlySeq_"+onlySeq + "CT "+cycTime);
		
		System.out.printf("%10s", "instance");
		System.out.printf("%10s", "run");
		System.out.printf("%10s", "LNS");
		System.out.printf("%14s", "time");
		System.out.printf("%14s", "z_FCFS");
		System.out.printf("%14s", "Gap-FCFS-LNS");
		System.out.println();
		if (csvProcess.getNumLines(csvFilePath)==0){
            ArrayList<String> header=new ArrayList<>();
            header.add("instance");
            header.add("run");
            header.add("LNS");
            header.add("time");
            header.add("z_FCFS");
            header.add("Gap-FCFS-LNS");
            csvProcess.writeHead(header,csvFilePath);
        }
		for (int i = start_ins; i < end_ins; i++) {
			Instance para = new Instance(numProd, numMach, numVeh, charge, vehSpeed, lineSpeed, cycTime,
			dataFile);
			para.Initialization(i);
			// 用fcfs 规则
			LogiTaskSortByTss logiTaskSortByTss = new LogiTaskSortByTss();
			ProdSortByWeightProcessTimeRatio prodSortByWeightProcessTimeRatio = new ProdSortByWeightProcessTimeRatio();
			ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
			for (int p = 0; p < para.numProd; p++) {
				ProdRatio prodRatio = new ProdRatio(p, 1 / ((para.numMach - 1 + para.demands[p]) * para.cycTime[p]
						+ para.moldChangeTime
						+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed));
				prodList.add(prodRatio);
			}
			Collections.sort(prodList, prodSortByWeightProcessTimeRatio);
			int[] prodSeq_fcfs = new int[para.numProd];
			for (int p = 0; p < para.numProd; p++) {
				prodSeq_fcfs[p] = prodList.get(p).getId();
			}
			ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
			String key;
			double timeCur = 0;// earliestReadyTime;// 时间游标
			for (int p = 0; p < para.numProd; p++) {// product从0开始编号
				for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
					if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
						timeCur = timeCur
								- para.prodSet.get(prodSeq_fcfs[p]).getCycleTime() * (para.demands[prodSeq_fcfs[p]] - 1)
								+ para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
					}
					for (int k = 0; k < para.K[prodSeq_fcfs[p]][m]; k++) {// batch从0开始编号
						key = prodSeq_fcfs[p] + "_" + m + "_" + k;
						// 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
						TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
						tssWithoutLogiConsList.add(tssWithoutLogiCons);
	//								PMKtoLogiTask.get(key).setTss(timeCur);
						timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
					}
				}
				timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq_fcfs[p]] + para.moldChangeTime
						- para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
			}
			Collections.sort(tssWithoutLogiConsList, logiTaskSortByTss);
			
			String[] pmk = new String[para.logiTaskSet.size()];
			for (int s = 0; s < para.logiTaskSet.size(); s++) {
				pmk[s] = tssWithoutLogiConsList.get(s).pmk;
			}
			ArrayList<ArrayList<Integer>> logiSeq_fcfs = new ArrayList<ArrayList<Integer>>();
			double[] timeCur2 = new double[para.numVeh];
			for (int r = 0; r < para.numVeh; r++) {
				logiSeq_fcfs.add(new ArrayList<Integer>());
			}
			for (int s = 0; s < para.logiTaskSet.size(); s++) {// for (int p = 0; p < prodSeq2.length; p++)
				// {//
				int logiTaskId = para.PMKtoLogiTask.get(pmk[s]).getId();// arrayOfI.get(index);//
				// 对应logiTask的ID
				int min_index = 0;
				double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
						+ para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
				double min_val = Float.MAX_VALUE;
				for (int r = 0; r < para.numVeh; r++) {
					if (timeCur2[r] + requiredDeliTime < min_val) {
						min_val = timeCur2[r] + requiredDeliTime;
						min_index = r;
					}
				}
				logiSeq_fcfs.get(min_index).add(logiTaskId);
				timeCur2[min_index] = timeCur2[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
			}
			// initialSol.printSol();
			// System.out.println("--------------");
			// for (int v = 0; v < logiSeq_fcfs.size(); v++) {
			// 	for (int j = 0; j < logiSeq_fcfs.get(v).size(); j++) {
			// 		System.out.print(logiSeq_fcfs.get(v).get(j));
			// 		System.out.print("->");
			// 	}
			// 	System.out.println();
			// }
			double z_FCFS=Float.MAX_VALUE;
			double bestObj=Float.MAX_VALUE;
			double bestbestUB = Float.MAX_VALUE;
			double st=0;
			double et=0;
			double time_meta=0;
			int jstar;
			for (int j = 0; j < 5; j++) {
				MetaHeur metaHeur = new MetaHeur(para,1,1,prodSeq_fcfs,logiSeq_fcfs,onlySeq);
				st = System.currentTimeMillis();
				
				bestObj=metaHeur.run(numProd,numVeh,para);
				
				et = System.currentTimeMillis();
				time_meta=(et-st)/1000;
				
				if (bestObj < bestbestUB) {
					bestbestUB = bestObj;
					jstar = j;
				}
				
				z_FCFS=metaHeur.z_FCFS;
				
				System.out.printf("%10d", i);
				System.out.printf("%10d", j);
				System.out.printf("%10.2f", bestObj);
				System.out.printf("%14.2f", time_meta);
				System.out.printf("%14.2f", z_FCFS);
				System.out.printf("%14.2f", 100*(z_FCFS-bestObj)/z_FCFS);
				System.out.println();
				ArrayList<String> line=new ArrayList<>();
			    line.add(String.valueOf(i));
			    line.add(String.valueOf(j));
			    line.add(String.valueOf(bestObj));
			    line.add(String.valueOf((et-st)/1000));
			    line.add(String.valueOf(z_FCFS));
			    line.add(String.valueOf(100*(z_FCFS-bestObj)/z_FCFS));
			    if (csvProcess.getNumLines(csvFilePath)==i*5+j+1) {
			    	csvProcess.writeALine(line, csvFilePath);
			    }else {
			    	csvProcess.replaceALine(line,csvFilePath,i*5+j+2);
				}
			}
		}

	}


}
