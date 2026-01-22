package subgradient_method;
import java.util.Comparator;

public class LogiTaskSortByPrior implements Comparator<LogiTaskPrior>{
    @Override
    public int compare(LogiTaskPrior o1, LogiTaskPrior o2) {
        if (o1.pri < o2.pri){
            return -1;
        } else if (o1.pri > o2.pri) {
            return 1;
        }
        else {
            return 0;
        }
    }
}
