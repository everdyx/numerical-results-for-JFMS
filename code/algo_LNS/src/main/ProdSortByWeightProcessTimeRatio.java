package main;

import java.util.Comparator;

public class ProdSortByWeightProcessTimeRatio implements Comparator<ProdRatio> {

	@Override
	public int compare(ProdRatio o1, ProdRatio o2) {
		return o1.getRatio() < o2.getRatio() ? 1 : -1;
	}

}
