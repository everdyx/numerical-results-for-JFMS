package basics;

import java.util.ArrayList;

public class ProdSolution {
	Instance para;
	double[][] xVal;
	ArrayList<Integer> prodSeq = new ArrayList<>();

	public ProdSolution(Instance para, double[][] xVal) {
		this.para = para;
		this.xVal = xVal;
	}


	public void parse() {
		// 模型可解，生成生产顺序
		boolean flag = true;
		int cur = para.numProd;
		while (flag) {
			for (int p = 0; p < para.numProd + 1; p++) {
				if (xVal[cur][p] >= 0.5) {
					prodSeq.add(p);
					cur = p;
					break;
				}
			}
			if (cur == para.numProd) {
				flag = false;
			}
		}
	}

	public void solCheck() {
		System.out.println("subtour check#########");
		if ( prodSeq.size() < para.numProd+1 ) {
			System.out.println("wrong solution");
		}
	}

	public void print() {
		System.out.println("production sequence:");
		for (int i=0;i<para.numProd;i++) {
			System.out.print(" -> "+prodSeq.get(i));
		}
		System.out.println(" ");
	}
}