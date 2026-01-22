package main;

import java.util.ArrayList;

import basics.Instance;

public class Solution {
	Instance para;
	double[][] xVal;
	ArrayList<Integer> prodSeq = new ArrayList<>();
	
	double[][] yVal;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();
	double[] tauVal;
	double[][][] tssVal;
	double[] tafVal;

	public Solution(Instance para, double[][] xVal, double[][] yVal, double[] tauVal, double[][][] tssVal, double[] tafVal) {
		this.para = para;
		this.xVal = xVal;
		this.yVal = yVal;
		this.tauVal=tauVal;
		this.tafVal=tafVal;
		this.tssVal=tssVal;
	}

	public void parse() {
		if (xVal!=null) {
		boolean flag1 = true;
		int cur = para.numProd;
		while (flag1) {
			for (int p = 0; p < para.numProd + 1; p++) {
				if (xVal[cur][p] >= 0.5) {
					prodSeq.add(p);
					cur = p;
					break;
				}
			}
			if (cur == para.numProd) {
				flag1 = false;
			}
		}
		}
		
		
		for (int k = 0; k < para.numVeh; k++) {
			ArrayList<Integer> r = new ArrayList<>(); // 定义一个对象为int型的链表
			routes.add(r); // 将上述定义的链表加入到链表routes中
		}
		// 模型可解，生成车辆路径
		for (int k = 0; k < para.numVeh; k++) {
			routes.get(k).add(-1);
		}
		int cnt = 0;
		for (int i = 0; i < yVal[para.logiTaskSet.size()].length; i++) {
			if (yVal[para.logiTaskSet.size()][i] >= 0.5) {
				routes.get(cnt).add(i);
				cnt++;
			}
		}
		for (int k = 0; k < para.numVeh; k++) {
			boolean flag2 = true;
			if (routes.get(k).size() <= 1) {
				continue;
			} else {
				int i = routes.get(k).get(1);
				while (flag2) {
					for (int j = 0; j < yVal[0].length; j++) {
						if (yVal[i][j] >= 0.5) {
							if (j==para.logiTaskSet.size()) {
								routes.get(k).add(-1);
							}
							else{
								routes.get(k).add(j);
							}
							i = j;
							break;
						}
					}
					if (i == para.logiTaskSet.size()) {
						flag2 = false;
					}
				}
			}
		}
	}

	public void print() {
		System.out.println("production sequence (Cplex):");
		for (int i=0;i<prodSeq.size();i++) {
			System.out.print(" -> "+prodSeq.get(i));
		}
		System.out.println();
		
		
		System.out.printf("%-16s\t", "prodID");
		for (int i = 0; i < para.prodSet.size(); i++) {
			System.out.printf("%10d", i);
		}
		System.out.println();
		
		System.out.printf("%-16s\t", "taf");
		for (int i = 0; i < para.prodSet.size(); i++) {
			System.out.printf("%10.4f", tafVal[i]);
		}
		System.out.println();
		
		//###############
		System.out.println("logistics sequence: ");
		for (int i=0;i<routes.size();i++) {
			for (int j=1;j<routes.get(i).size()-1;j++) {
				if (routes.get(i).get(j)<para.logiTaskSet.size()) {
					System.out.print(" -> "+routes.get(i).get(j)+"("+tauVal[routes.get(i).get(j)]+",P"+para.ItoLogiTask.get(routes.get(i).get(j)).getProd().getId()+",M"+para.ItoLogiTask.get(routes.get(i).get(j)).getMach().getId()+")");
				}
			}
			System.out.println();
		}
//		for (int i=0;i<routes.size();i++) {
//			for (int j=1;j<routes.get(i).size()-1;j++) {
//				if(routes.get(i).get(j)<para.logiTaskSet.size()) {
//					System.out.print(" -> M"+para.ItoLogiTask.get(routes.get(i).get(j)).getMach().getId());
//				}else {
//					System.out.print(" -> Charing");
//				}
//			}
//			System.out.println();
//		}
		
		System.out.printf("%-16s\t", "LogiID");
		for (int i = 0; i < para.logiTaskSet.size(); i++) {
			System.out.printf("%10d", i);
		}
		System.out.println();

		System.out.printf("%-16s\t", "tad");
		for (int i = 0; i < para.logiTaskSet.size(); i++) {
			System.out.printf("%10.4f", tauVal[i]);
		}
		System.out.println();
		
		System.out.printf("%-16s\t", "tss");
		double[] temp=new double[para.logiTaskSet.size()];
		for (int p = 0; p < para.prodSet.size(); p++) {
			for (int m = 0; m < para.numMach; m++) {
				for (int k = 0; k < para.K[p][m]; k++) {
					temp[para.PMKtoLogiTask.get(p+"_"+m+"_"+k).getId()]=tssVal[p][m][k];
					System.out.printf("%10.4f", temp[para.PMKtoLogiTask.get(p+"_"+m+"_"+k).getId()]);
				}
			}
		}
//		for (int i = 0; i < para.logiTaskSet.size(); i++) {
//			System.out.printf("%10.4f", temp[i]);
//		}
		System.out.println();
//		System.out.println(routes);
	}

}
