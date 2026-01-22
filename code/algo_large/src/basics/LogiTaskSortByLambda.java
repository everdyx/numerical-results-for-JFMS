package basics;

import java.util.Comparator;

public class LogiTaskSortByLambda implements Comparator<LogiTaskIdLambda> {
    @Override
    public int compare(LogiTaskIdLambda o1, LogiTaskIdLambda o2) {
        if(o1.lambda<o2.lambda){
            return -1;
        } else if (o1.lambda==o2.lambda) {
            return 0;
        }else{
            return 1;
        }
    }
}
