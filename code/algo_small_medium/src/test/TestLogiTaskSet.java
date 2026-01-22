package test;

import java.util.Arrays;
import java.util.HashMap;

public class TestLogiTaskSet {
	public static void main(String[] args) {
//		Instance ins=new Instance();
//		ins.Initialization();
//		String key = 0 + "_" + 0 + "_" + 0;
//		System.out.println(ins.logiTaskSet);
//		System.out.println(ins.logiTaskSet.get(0).getProd().getId());
//		System.out.println(ins.PMKtoLogiTask.get(key));
//		Collections.sort(ins.logiTaskSet, new LogiTaskSortByLambdaTransTimeRatio());
//		System.out.println(ins.logiTaskSet);
//		System.out.println(ins.logiTaskSet.get(0).getProd().getId());
//		System.out.println(ins.PMKtoLogiTask.get(key));
//		TestTableOutput a=new TestTableOutput();
//		System.out.println(a.a(0));
		int[] prodSeq= {1,2,0};
		int[] prodSeq2= {1,2,0};
		String str1=Arrays.toString(prodSeq);
		String str2=Arrays.toString(prodSeq);
		HashMap<String, Double> mapHeur3=new HashMap<String, Double>();
		mapHeur3.put(str1, 3.0);
		mapHeur3.put(str2, 6.0);
		System.out.println(mapHeur3.get(str1));
	}
	

}
