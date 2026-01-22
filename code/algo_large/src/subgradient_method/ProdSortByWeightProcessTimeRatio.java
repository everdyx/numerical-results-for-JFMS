package subgradient_method;

import java.util.Comparator;

public class ProdSortByWeightProcessTimeRatio implements Comparator<ProdRatio> {

	@Override
	public int compare(ProdRatio o1, ProdRatio o2) {
		if(o1.getRatio() < o2.getRatio()){
			return 1;
		}else if(o1.getRatio() == o2.getRatio()){
			return 0;
		}else {
			return -1;
		}
	}

}
