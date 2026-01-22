package basics;

import java.util.ArrayList;

public class LogiSolution_SumC {
	Instance para;
	
	double[][] yVal;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();
	double[] tauVal;


	public LogiSolution_SumC(Instance para, double[][] yVal, double[] tauVal) {
		this.para = para;
		this.yVal = yVal;
		this.tauVal=tauVal;
	}

	public void parse() {
		
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
		
		//###############
		System.out.println("logistics sequence: ");
		for (int i=0;i<routes.size();i++) {
			for (int j=1;j<routes.get(i).size()-1;j++) {
				System.out.print(" -> "+routes.get(i).get(j));
			}
			System.out.println();
		}
		for (int i=0;i<routes.size();i++) {
			for (int j=1;j<routes.get(i).size()-1;j++) {
				if(routes.get(i).get(j)<para.logiTaskSet.size()) {
					System.out.print(" -> M"+para.ItoLogiTask.get(routes.get(i).get(j)).getMach().getId());
				}else {
					System.out.print(" -> Charing");
				}
			}
			System.out.println();
		}
		
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
		
//		System.out.println(routes);
	}

}
