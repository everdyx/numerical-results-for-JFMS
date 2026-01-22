package basics;

import java.util.ArrayList;

public class LogiSolution {
	Instance para;
	double[][] yVal;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();

	public LogiSolution(Instance para, double[][] yVal) {
		this.para = para;
		this.yVal = yVal;
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
//		if (cnt != para.numVeh) {
//			System.out.println("only " + cnt + "vehicles are used, wrong solution!!!");
//		}
		for (int k = 0; k < para.numVeh; k++) {
			boolean flag = true;
			if (routes.get(k).size() <= 1) {
				continue;
			} else {
				int i = routes.get(k).get(1);
				while (flag) {
					for (int j = 0; j < para.logiTaskSet.size()+1; j++) {
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
						flag = false;
					}
				}
			}
		}
	}

	public void solCheck() {
		// TODO Auto-generated method stub
		System.out.println("subtour check#########");
		ArrayList<Integer> visitedVertex = new ArrayList<>();
		ArrayList<Integer> notVisitedVertex = new ArrayList<>();
		for (int k = 0; k < routes.size(); k++) {
			for (int i = 1; i < routes.get(k).size() - 1; i++) {
				visitedVertex.add(routes.get(k).get(i));
			}
		}
		for (int i = 0; i < para.logiTaskSet.size(); i++) {
			if (!visitedVertex.contains(i)) {
				notVisitedVertex.add(i);
				for (int j = 0; j < yVal[i].length; j++) {
					if (yVal[i][j] > 0.5) {
						System.out.println(i + "in subtour");
					}
				}
			}
		}
	}
	
	public void print() {
		System.out.println(routes);
	}
}
