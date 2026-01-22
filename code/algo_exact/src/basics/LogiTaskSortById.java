package basics;

import java.util.Comparator;

public class LogiTaskSortById implements Comparator<LogiTask>{

	@Override
	public int compare(LogiTask o1, LogiTask o2) {
		return o1.getId()>o2.getId()?1:-1;
	}

}
