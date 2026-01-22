package basics;

import java.util.Comparator;

public class LogiTaskSortByLambdaProcessTimeRatio implements Comparator<LogiTask>{

    @Override
    public int compare(LogiTask o1, LogiTask o2) {
        if ((o1.getLambda()/o1.getProcessTime())<(o2.getLambda()/o2.getProcessTime())){
            return 1;
        } else if ((o1.getLambda()/o1.getProcessTime())==(o2.getLambda()/o2.getProcessTime())) {
            return 0;
        }else {
            return -1;
        }
    }

}
