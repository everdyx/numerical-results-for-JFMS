package basics;

import java.util.Comparator;

public class LogiTaskSortByLambdaTransTimeRatio implements Comparator<LogiTask>{

	@Override
	public int compare(LogiTask o1, LogiTask o2) {
		return (o1.getLambda()/o1.getTransTime())<(o2.getLambda()/o2.getTransTime())?1:-1;
	}

}
