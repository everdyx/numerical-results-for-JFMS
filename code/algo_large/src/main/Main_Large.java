package main;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;

import basics.Instance;
import subgradient_method.LogiTaskSortByTss;
import subgradient_method.ProdRatio;
import subgradient_method.ProdSortByWeightProcessTimeRatio;
import subgradient_method.SubgradMethod;
import subgradient_method.TssWithoutLogiCons;

public class Main_Large {

	public void run(int start_ins, int end_ins, int numProd, int numVeh) throws Exception {
//		int start_ins = 0;
//		int end_ins = start_ins+100;
		boolean charge;
//		int numProd = 10;
		double cycTime = 3;
		double vehSpeed = 1;// m/s
		boolean newLB = false;
		double stepsizeFactor = 1;// 一般在0到2之间

		if (numProd <= 3) {
			charge = false;
		} else {
			charge = true;
		}

		int numMach = 8;
//		int numVeh = 2 ;

		String csvFilePath = ".\\output\\" + "P" + numProd + "M" + numMach + "V" + numVeh + " CT" + cycTime + ".csv";
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
		System.out.println("P" + numProd + "M" + numMach + "V" + numVeh + " CT" + cycTime);

		double lineSpeed = 0.5;
		System.out.println("P" + numProd + "M" + numMach + "V" + numVeh);

		if (csvProcess.getNumLines(csvFilePath) == 0) {
			ArrayList<String> header = new ArrayList<>();
			header.add("instance");
			header.add("run");
			header.add("LB");
			header.add("UB");
			header.add("FCFS");
			header.add("time");
			header.add("G(UB-LB)");
			header.add("G(FCFS-UB)");
			header.add("LBNoLift");
			header.add("UBNoLS");
			header.add("LBImprov");
			header.add("UBImprov");
			header.add("prodSeq");
			header.add("logiSeq");
			csvProcess.writeHead(header, csvFilePath);
		}

		Instance para = new Instance(numProd, numMach, numVeh, charge, vehSpeed, lineSpeed, cycTime,
				"data\\demand3.txt");
		ProdSortByWeightProcessTimeRatio prodSortByWeightProcessTimeRatio = new ProdSortByWeightProcessTimeRatio();
		LogiTaskSortByTss logiTaskSortByTss = new LogiTaskSortByTss();
		System.out.printf("%10s", "instance");
		System.out.printf("%10s", "run");
		System.out.printf("%10s", "LB");
		System.out.printf("%10s", "UB");
		System.out.printf("%10s", "FCFS");
		System.out.printf("%14s", "time");
		System.out.printf("%14s", "G(UB-LB)");
		System.out.printf("%14s", "G(FCFS-UB)");
		System.out.printf("%14s", "LBNoLift");
		System.out.printf("%14s", "UBNoLS");
		System.out.printf("%14s", "LBImprov");
		System.out.printf("%14s", "UBImprov");
		System.out.println();
		for (int i = start_ins; i < end_ins; i++) {
			para.Initialization(i);

			// FCFS
			SubgradMethod subgradMethodFCFS = new SubgradMethod(para, stepsizeFactor, newLB);
			subgradMethodFCFS.initialization();
			// 先生成生产顺序
			ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
			for (int p = 0; p < para.numProd; p++) {
				ProdRatio prodRatio = new ProdRatio(p, 1 / ((para.numMach - 1 + para.demands[p]) * para.cycTime[p]
						+ para.moldChangeTime
						+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed));
				prodList.add(prodRatio);
			}
			Collections.sort(prodList, prodSortByWeightProcessTimeRatio);
			int[] prodSeq2 = new int[para.numProd];
			for (int p = 0; p < para.numProd; p++) {
				prodSeq2[p] = prodList.get(p).getId();
			}
			// 后物流顺序
			// 如果是先呼叫先服务策略，物流顺序由tss的先后次序决定
			// 先得到tss,并记录对应的pmk
			ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
			String key;
			double timeCur = 0;// earliestReadyTime;// 时间游标
			for (int p = 0; p < para.numProd; p++) {// product从0开始编号
				for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
					if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
						timeCur = timeCur
								- para.prodSet.get(prodSeq2[p]).getCycleTime() * (para.demands[prodSeq2[p]] - 1)
								+ para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
					}
					for (int k = 0; k < para.K[prodSeq2[p]][m]; k++) {// batch从0开始编号
						key = prodSeq2[p] + "_" + m + "_" + k;
						// 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
						TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
						tssWithoutLogiConsList.add(tssWithoutLogiCons);
//								PMKtoLogiTask.get(key).setTss(timeCur);
						timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
					}
				}
				timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq2[p]] + para.moldChangeTime
						- para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
			}
			Collections.sort(tssWithoutLogiConsList, logiTaskSortByTss);
			String[] pmk = new String[para.logiTaskSet.size()];
			for (int s = 0; s < para.logiTaskSet.size(); s++) {
				pmk[s] = tssWithoutLogiConsList.get(s).pmk;
			}
			ArrayList<ArrayList<Integer>> logiSeq = new ArrayList<ArrayList<Integer>>();
			double[] timeCur2 = new double[para.numVeh];
			for (int r = 0; r < para.numVeh; r++) {
				logiSeq.add(new ArrayList<Integer>());
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
				logiSeq.get(min_index).add(logiTaskId);
				timeCur2[min_index] = timeCur2[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
			}
			double z_FCFS = subgradMethodFCFS.getObj(prodSeq2, logiSeq);
			int[] prodSeq_f=prodSeq2.clone();
			ArrayList<ArrayList<Integer>> logiSeq_f=subgradMethodFCFS.listCopy2(logiSeq);

			double st, et;
			double subgradMethodTime = 0;
			double bestUB = Float.MAX_VALUE;
			double bestbestUB = Float.MAX_VALUE;
			double bestUBNoLS = Float.MAX_VALUE;
			double bestLB = Float.MIN_VALUE;
			double bestLBNoLift = Float.MIN_VALUE;
			int jstar = 0;
			int[][] prodSeq_best = new int[5][para.numProd];
			ArrayList<ArrayList<ArrayList<Integer>>> logiSeq_best = new ArrayList<ArrayList<ArrayList<Integer>>>();
			// 次梯度算法
			for (int j = 0; j < 5; j++) {
				SubgradMethod subgradMethod = new SubgradMethod(para, stepsizeFactor, newLB);
				subgradMethod.initialization();
				st = System.currentTimeMillis();
				
				subgradMethod.run_surrogate(z_FCFS, prodSeq2, logiSeq);
				
				et = System.currentTimeMillis();
				subgradMethodTime = (et - st) / 1000;
				bestLBNoLift = subgradMethod.bestRelaxObj;
				bestLB = subgradMethod.LBLift;
				
				bestUB = subgradMethod.getObj(subgradMethod.bestProdSeq, subgradMethod.bestLogiSeq);
//				double aa=subgradMethodFCFS.getObj(prodSeq2, logiSeq);
//				double bb=subgradMethod.getObj(prodSeq2, logiSeq);
//				double cc=subgradMethod.getObj(prodSeq_f, logiSeq_f);
				if(bestUB>z_FCFS) {
					bestUB=z_FCFS;
					prodSeq_best[j]=prodSeq_f.clone();
					logiSeq_best.add(logiSeq_f);
				}else {
					prodSeq_best[j] = subgradMethod.bestProdSeq.clone();
					logiSeq_best.add(subgradMethod.bestLogiSeq);
				}
				
				if (bestUB< bestbestUB) {
					bestbestUB = bestUB;
					jstar = j;
				}
				bestUBNoLS = subgradMethod.UBNoImprove;

				System.out.printf("%10d", i);
				System.out.printf("%10d", j);
				System.out.printf("%10.2f", bestLB);
				System.out.printf("%10.2f", bestUB);
				System.out.printf("%10.2f", z_FCFS);
				System.out.printf("%14.2f", subgradMethodTime);
				System.out.printf("%14.2f", 100 * (bestUB - bestLB) / bestUB);
				System.out.printf("%14.2f", 100 * (z_FCFS - bestUB) / z_FCFS);
				System.out.printf("%14.2f", bestLBNoLift);// System.out.printf("%14s", "LBNoLift");
				System.out.printf("%14.2f", bestUBNoLS);// System.out.printf("%14s", "UBNoLS");
				System.out.printf("%14.2f", 100 * (bestLB - bestLBNoLift) / bestLB);// System.out.printf("%14s",
																					// "LBImprov");
				System.out.printf("%14.2f", 100 * (bestUBNoLS - bestUB) / bestUBNoLS);// System.out.printf("%14s",
				System.out.printf("%14.2f", subgradMethod.getObj(prodSeq_best[jstar], logiSeq_best.get(jstar)));
				System.out.printf("%14.2f",
						subgradMethod.getObj(prodSeq_best[jstar], logiSeq_best.get(jstar)) - bestUB);																		// "UBImprov");
				System.out.println();
				if( Math.abs( subgradMethod.getObj(prodSeq_best[j], logiSeq_best.get(j)) - bestUB)>1 ) {
					System.out.println("wrong");
					System.out.println("subg"+bestUB);
					System.out.println("best"+subgradMethod.getObj(prodSeq_best[j], logiSeq_best.get(j)));
					System.exit(1);
				}
				
				ArrayList<String> line = new ArrayList<>();
				line.add(String.valueOf( i));
				line.add(String.valueOf(j));
				line.add(String.valueOf(bestLB));
				line.add(String.valueOf(bestUB));
				line.add(String.valueOf(z_FCFS));
				line.add(String.valueOf(subgradMethodTime));
				line.add(String.valueOf(100 * (bestUB - bestLB) / bestUB));
				line.add(String.valueOf(100 * (z_FCFS - bestUB) / z_FCFS));
				line.add(String.valueOf(bestLBNoLift));
				line.add(String.valueOf(bestUBNoLS));
				line.add(String.valueOf(100 * (bestLB - bestLBNoLift) / bestLB));
				line.add(String.valueOf(100 * (bestUBNoLS - bestUB) / bestUBNoLS));
				
				int[] array = prodSeq_best[jstar];
				String result1 = Arrays.stream(array).mapToObj(String::valueOf).collect(Collectors.joining("-"));
				line.add(result1);

				ArrayList<ArrayList<Integer>> listOfLists = logiSeq_best.get(jstar);
				String result2 = listOfLists.stream()
						.map(innerList -> innerList.stream().map(String::valueOf).collect(Collectors.joining("-")))
						.collect(Collectors.joining(" | "));
				line.add(result2);
				
				if (csvProcess.getNumLines(csvFilePath) == i * 5 + j + 1) {
					csvProcess.writeALine(line, csvFilePath);
				} else {
					csvProcess.replaceALine(line, csvFilePath, i * 5 + j + 2);
				}
			}
		}
	}

}
