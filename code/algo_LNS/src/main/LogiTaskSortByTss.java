package main;

import java.util.Comparator;

public class LogiTaskSortByTss implements Comparator<TssWithoutLogiCons>{
	@Override
	public int compare(TssWithoutLogiCons o1, TssWithoutLogiCons o2) {
		return o1.tss < o2.tss ? -1 : 1;
	}
}