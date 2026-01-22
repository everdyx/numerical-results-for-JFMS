package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Scanner;
import java.util.stream.Collectors;

import basics.Instance;
import ilog.cplex.IloCplex;
import subgradient_method.LogiTaskSortByTss;
import subgradient_method.ProdRatio;
import subgradient_method.ProdSortByWeightProcessTimeRatio;
import subgradient_method.SubgradMethod;
import subgradient_method.TssWithoutLogiCons;

public class Main_Exact {

	public void run(int start_ins, int end_ins, int numProd, int numMach, int numVeh, double cycTime, boolean charge,
			boolean given_q, double timeLimit) throws Exception {
//        int start_ins =0 ;
//        int end_ins = start_ins+100;
//        boolean charge=false;
//        int numProd = 2;
//        double cycTime = 3;
		String csvFilePath = ".\\output\\" + "P" + numProd + "M" + numMach + "V" + numVeh + " CT" + cycTime + "_charge_"
				+ charge + "_giveQ_" + given_q + "_time_limit_" + timeLimit + ".csv";
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
		CSVProcess csvProcess = new CSVProcess();
		System.out.println("P" + numProd + "M" + numMach + "V" + numVeh + " CT" + cycTime + "_charge_" + charge
				+ "_giveQ_" + given_q);
		double vehSpeed = 1;// m/s

//        int numMach = 2;
//        int numVeh = 2;
//        boolean given_q = false;
		double lineSpeed = 0.5;

		ProdSortByWeightProcessTimeRatio prodSortByWeightProcessTimeRatio = new ProdSortByWeightProcessTimeRatio();
		System.out.printf("%10s", "instance");
		System.out.printf("%14s", "SD");
		System.out.printf("%14s", "SDgap");
		System.out.printf("%14s", "E");
		System.out.printf("%14s", "Egap");
		System.out.printf("%14s", "time_SD");
		System.out.printf("%14s", "time_E");
		System.out.printf("%14s", "Gap-SD-E");
		System.out.println();
		if (csvProcess.getNumLines(csvFilePath) == 0) {
			ArrayList<String> header = new ArrayList<>();
			header.add("instance");
			header.add("SD");
			header.add("SDgap");
			header.add("E");
			header.add("Egap");
			header.add("time_SD");
			header.add("time_E");
			header.add("Gap-SD-E");
			header.add("prodSeq");
			header.add("logiSeq");
			csvProcess.writeHead(header, csvFilePath);
		}
		Instance para;
		LogiTaskSortByTss logiTaskSortByTss = new LogiTaskSortByTss();

		if (numMach > 2) {
			para = new Instance(numProd, numMach, numVeh, charge, vehSpeed, lineSpeed, cycTime, "data\\demand3.txt");
		} else {
			para = new Instance(numProd, numMach, numVeh, charge, vehSpeed, lineSpeed, cycTime, "data\\demand2.txt");
		}
		for (int i = start_ins; i < end_ins; i++) {
			para.Initialization(i);
			double st, et;
			// 启发式
			SubgradMethod subgradMethod_fcff = new SubgradMethod(para, 1, false);
			subgradMethod_fcff.initialization();
			st = System.currentTimeMillis();
			// 先生成生产顺序
			ArrayList<ProdRatio> prodList_f = new ArrayList<ProdRatio>();
			for (int p = 0; p < para.numProd; p++) {
				ProdRatio prodRatio = new ProdRatio(p, 1 / ((para.numMach - 1 + para.demands[p]) * para.cycTime[p]
						+ para.moldChangeTime
						+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed));
				prodList_f.add(prodRatio);
			}
			Collections.sort(prodList_f, prodSortByWeightProcessTimeRatio);
			int[] prodSeq_fcff = new int[para.numProd];
			for (int p = 0; p < para.numProd; p++) {
				prodSeq_fcff[p] = prodList_f.get(p).getId();
			}
			// 物流顺序采用先呼叫先服务策略，物流顺序由tss的先后次序决定
			// 先得到tss,并记录对应的pmk
			ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
			String key;
			double timeCur = 0;// earliestReadyTime;// 时间游标
			for (int p = 0; p < para.numProd; p++) {// product从0开始编号
				for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
					if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
						timeCur = timeCur
								- para.prodSet.get(prodSeq_fcff[p]).getCycleTime() * (para.demands[prodSeq_fcff[p]] - 1)
								+ para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
					}
					for (int k = 0; k < para.K[prodSeq_fcff[p]][m]; k++) {// batch从0开始编号
						key = prodSeq_fcff[p] + "_" + m + "_" + k;
						// 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
						TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
						tssWithoutLogiConsList.add(tssWithoutLogiCons);
						timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
					}
				}
				timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq_fcff[p]] + para.moldChangeTime
						- para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
			}
			Collections.sort(tssWithoutLogiConsList, logiTaskSortByTss);
			String[] pmk = new String[para.logiTaskSet.size()];
			for (int s = 0; s < para.logiTaskSet.size(); s++) {
				pmk[s] = tssWithoutLogiConsList.get(s).pmk;
			}
			ArrayList<ArrayList<Integer>> logiSeq_fcff = new ArrayList<ArrayList<Integer>>();
			double[] timeCur2 = new double[para.numVeh];
			for (int r = 0; r < para.numVeh; r++) {
				logiSeq_fcff.add(new ArrayList<Integer>());
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
				logiSeq_fcff.get(min_index).add(logiTaskId);
				timeCur2[min_index] = timeCur2[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
			}
			double obj_fcff = subgradMethod_fcff.getObj(prodSeq_fcff, logiSeq_fcff);
			et = System.currentTimeMillis();
			double time_fcff = (et - st) / 1000;

//			求解器求解
			double exactTime1 = 0;
			double exactTime = 0;
			PrimalProblem primalProblem = new PrimalProblem(para, given_q, charge,obj_fcff);
			primalProblem.show = false;
			ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
			for (int p = 0; p < para.numProd; p++) {
				ProdRatio prodRatio = new ProdRatio(p, 1 / ((para.numMach - 1 + para.demands[p]) * para.cycTime[p]));
				prodList.add(prodRatio);
			}
			Collections.sort(prodList, prodSortByWeightProcessTimeRatio);
			int[] prodSeq2 = new int[para.numProd];
			for (int p = 0; p < para.numProd; p++) {
				prodSeq2[p] = prodList.get(p).getId();
			}
			PrimalProblem primalProblem1 = new PrimalProblem(para, prodSeq2, given_q, charge,obj_fcff);
			primalProblem1.show = false;
//			System.out.println("===========solution without collaboration============");
			st = System.currentTimeMillis();
			primalProblem1.buildModel();
			primalProblem1.model.setParam(IloCplex.Param.TimeLimit, timeLimit);
			primalProblem1.getSolution();
			et = System.currentTimeMillis();
			exactTime1 = (et - st) / 1000;
//			System.out.println("===========solution with collaboration============");
			st = System.currentTimeMillis();
			primalProblem.buildModel();
			primalProblem.model.setParam(IloCplex.Param.TimeLimit, timeLimit);
			primalProblem.getSolution();
			et = System.currentTimeMillis();
			exactTime = (et - st) / 1000;

			System.out.printf("%10d", i);
			System.out.printf("%14.2f", primalProblem1.getObjVal());
			System.out.printf("%14.2f", primalProblem1.model.getMIPRelativeGap());
			System.out.printf("%14.2f", primalProblem.getObjVal());
			System.out.printf("%14.2f", primalProblem.model.getMIPRelativeGap());
			System.out.printf("%14.2f", exactTime1);
			System.out.printf("%14.2f", exactTime);
			System.out.printf("%14.2f",
					100 * (primalProblem1.getObjVal() - primalProblem.getObjVal()) / primalProblem1.getObjVal());
			System.out.println();

			ArrayList<String> line = new ArrayList<>();
			line.add(String.valueOf(i));
			line.add(String.valueOf(primalProblem1.getObjVal()));
			line.add(String.valueOf(primalProblem1.model.getMIPRelativeGap()));
			line.add(String.valueOf(primalProblem.getObjVal()));
			line.add(String.valueOf(primalProblem.model.getMIPRelativeGap()));
			line.add(String.valueOf(exactTime1));
			line.add(String.valueOf(exactTime));
			line.add(String.valueOf(
					100 * (primalProblem1.getObjVal() - primalProblem.getObjVal()) / primalProblem1.getObjVal()));

			primalProblem.getSolution();
			ArrayList<Integer> prodSeq = primalProblem.sol.prodSeq;
			int[] prodSeq_array = new int[para.numProd];
			for (int j = 0; j < prodSeq_array.length; j++) {
				prodSeq_array[j] = prodSeq.get(j);
			}
			
			String result1 = Arrays.stream(prodSeq_array).mapToObj(String::valueOf).collect(Collectors.joining("-"));
			line.add(result1);

			ArrayList<ArrayList<Integer>> logiSeq = primalProblem.sol.routes;
			for (int j = 0; j < logiSeq.size(); j++) {
				logiSeq.get(j).remove(logiSeq.get(j).size()-1);
				logiSeq.get(j).remove(0);
			}
			String result2 = logiSeq.stream()
					.map(innerList -> innerList.stream().map(String::valueOf).collect(Collectors.joining("-")))
					.collect(Collectors.joining(" | "));
			line.add(result2);
			
			if( Math.abs( subgradMethod_fcff.getObj(prodSeq_array, logiSeq) - primalProblem.getObjVal())>1 ) {
				System.out.println("wrong");
				System.exit(1);
			}

			if (csvProcess.getNumLines(csvFilePath) == i + 1) {
				csvProcess.writeALine(line, csvFilePath);
			} else {
				csvProcess.replaceALine(line, csvFilePath, i + 2);
			}
		}

	}
}
