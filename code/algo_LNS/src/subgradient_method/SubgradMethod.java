package subgradient_method;

import java.io.IOException;
import java.util.*;

import basics.*;
import column_generation.ColGen;
import ilog.concert.IloException;
import main.PrimalProblem;

public class SubgradMethod {

//    public double[] lambda_sur;
//    public double[] tad_sur;
//    public double[] tss_sur;
//    public double[] taf_sur;

    public double UBNoImprove=Float.MAX_VALUE;
    public double LBLift =Float.MIN_VALUE;
    boolean show = false;
    public Instance para;
    public double numerator_for_stepsize = 0.0002;
    public int iter = 0;
    public boolean infoOutput = false;
    public int iter_for_stepsize_factor = 0;
    public int iter_no_improve = 0;
    public int max_iter_no_improve = 10;
    public int max_iter_no_improve_for_stepsize_factor = 5;
    public double stepsize = Float.MAX_VALUE;
    public double stepsizeFactor;
    double[] stepsizeUpperBound_Prod;
    public double[] subgrad;
    //	public double[] lambda;
//	public double[] bestLambdaForUpperBound;
    public double[] bestLambdaForLowerBound;
    public double bestRelaxObj = Float.MIN_VALUE;
    public double prodObj;
    public double logiObj;
    public double earliestReadyTime = Float.MAX_VALUE;// 所有产品共同的最早开始时间
    public double makespan = 0;// 所有产品共同的最早开始时间
    double precision = 1e-4; // 0.0001 设置算法终止的精度
    public double objUpBound = Float.MAX_VALUE;
    double[] tad;
    double[] taf;
    public double bestObjUpBound = Float.MAX_VALUE;
    public int bestIter;
    private boolean improved;
    public double LBNoLift;
    public ArrayList<AgvPath> optAgvPaths = new ArrayList<AgvPath>();
    public ArrayList<LogiTask> logiTaskSet = new ArrayList<LogiTask>();
    ArrayList<ArrayList<Integer>> logiSeq = new ArrayList<ArrayList<Integer>>();
    ArrayList<ArrayList<Integer>> logiSeqHeur3 = new ArrayList<ArrayList<Integer>>();
    //	public boolean compare=false;
    public int[] prodSeq;
    public int[] bestProdSeq;
    double[] tard;
    double minTard;
    double[] tardLB;
    double[] tard_j;
    double[] tard2;
    public HashMap<String, LogiTask> PMKtoLogiTask = new HashMap<String, LogiTask>();
    //	HashMap<String, Double> mapHeur2=new HashMap<String, Double>();
//	HashMap<String, ArrayList<ArrayList<Integer>> > mapHeurProdSeq_LogiSeq=new HashMap<String, ArrayList<ArrayList<Integer>> >();
    public ArrayList<ArrayList<Integer>> bestLogiSeq = new ArrayList<ArrayList<Integer>>();
    private boolean newLB;
    private double maxTaf;
    private double maxTard;
    private double[] tafUB;
    public int[] iniProdSeq;
    public double prob;
    private double[] tafLB;

    public SubgradMethod(Instance para, double stepsizeFactor, boolean newLB) {
        this.para = para;
        this.stepsizeFactor = stepsizeFactor;
        this.newLB = newLB;
    }

    public void initialization() throws IloException, ClassNotFoundException, IOException {
        subgrad = new double[para.numProd];
        prodSeq = new int[para.numProd];
        iniProdSeq = new int[para.numProd];
        taf = new double[para.numProd];

        tad = new double[para.logiTaskSet.size()];
        bestProdSeq = new int[para.numProd];
        this.PMKtoLogiTask = para.PMKtoLogiTask;
        this.logiTaskSet = para.logiTaskSet;
        Collections.sort(logiTaskSet, new LogiTaskSortById());

        stepsizeUpperBound_Prod = new double[para.numProd];// stepSize的上界，确保新的lambda求和不会超过1
        for (int p = 0; p < para.numProd; p++) {
            stepsizeUpperBound_Prod[p] = Float.MAX_VALUE;
        }

        // ###计算earliestReadyTime
        earliestReadyTime = 0;
        double[] denominator_wrt_prod = new double[para.numProd];

        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p, (para.objValFactor) / denominator_wrt_prod[p]);
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());
        int[] prodSeq = new int[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
            iniProdSeq[p] = prodList.get(p).getId();
//			System.out.print("-> " + prodSeq[p] + " ");
        }
        /// 计算一个基础时间/////////////////////////////////////////
        if(para.numMach==2&&para.numVeh==2) {
            tard = new double[para.numProd];
            minTard = Float.MAX_VALUE;
            for (int p = 0; p < para.numProd; p++) {
                PrimalProblem primalProblem1 = new PrimalProblem(para, prodSeq[p]);
                primalProblem1.buildModel_1Prod();
                double productionTimeWithoutLogiConstriants = para.prodSet.get(prodSeq[p]).getCycleTime()
                        * (para.numMach - 1 + para.demands[prodSeq[p]])
                        + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
                tard[prodSeq[p]] = primalProblem1.model.getObjValue() - productionTimeWithoutLogiConstriants;
                if (show) {
                    System.out.println(prodSeq[p] + "_" + tard[prodSeq[p]] + "=" + primalProblem1.model.getObjValue() + "-"
                            + productionTimeWithoutLogiConstriants);
                }
//                if (tard[prodSeq[p]] < minTard) {
//                    minTard = tard[prodSeq[p]];
//                }
            }
            minTard=Arrays.stream(tard).min().getAsDouble();
        }else {
            minTard=0;
        }
//		double minTimePerLogiTask = para.getDistance(para.depot, para.machSet.get(0)) * 2 * para.consumingRate
//				* para.chargingRate + para.getDistance(para.depot, para.machSet.get(0)) * 2 / para.vehSpeed;
//		double T_LO = 0;
//		tard_j = new double[para.numProd];
//		tard_j[0] = minTard;
//		double T_PR = tard_j[0]
//				+ para.prodSet.get(prodSeq[0]).getCycleTime() * (para.numMach - 1 + para.demands[prodSeq[0]])
//				+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
//		for (int p = 1; p < para.numProd; p++)
//		{
//			PrimalProblem primalProblem1 = new PrimalProblem(para, prodSeq[p], T_PR, T_LO);
//			primalProblem1.buildModel_1Prod_v2();
//			double productionTimeWithoutLogiConstriants = para.prodSet.get(prodSeq[p]).getCycleTime()
//					* (para.numMach - 1 + para.demands[prodSeq[p]])
//					+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
//			double tardxx = primalProblem1.model.getObjValue() - T_PR - productionTimeWithoutLogiConstriants;
////			if (tardxx>=0.0) {
//			T_LO = T_LO + Math.min(para.K[p][0] / para.numVeh, 1) * minTimePerLogiTask;
////			}
//			T_PR = primalProblem1.model.getObjValue();
//			tard_j[p] = tardxx;
//		}
///////////////////////////////////////////////////////////
        // 计算maxTaf
//		for (int p = 0; p < para.numProd; p++)
//		{
//			prodSeq[p] = prodList.get(para.numProd - 1 - p).getId();
//		}
        double timeCur = minTard;

        // 拿一个初始解，通过初始解拿到最大的taf
        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
        }
        ArrayList<ArrayList<Integer>> logiSeq = new ArrayList<ArrayList<Integer>>();
        double[] timeCurVeh = new double[para.numVeh];
        double[] tad = new double[para.logiTaskSet.size()];
        for (int r = 0; r < para.numVeh; r++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }

            for (int index = 0; index < arrayOfI.size(); index++) {
                int logiTaskId = arrayOfI.get(index);// 对应logiTask的ID
                int min_index = 0;
                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                double min_val = Float.MAX_VALUE;
                for (int r = 0; r < para.numVeh; r++) {
                    if (timeCurVeh[r] + requiredDeliTime < min_val) {
                        min_val = timeCurVeh[r] + requiredDeliTime;
                        min_index = r;
                    }
                }
                logiSeq.get(min_index).add(logiTaskId);
                tad[logiTaskId] = timeCurVeh[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCurVeh[min_index] = timeCurVeh[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }
        }
        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur;
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        timeCur = 0;// 时间游标
        double sum_taf = 0;
        double[] taf = new double[para.numProd];
        key = prodSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(para.moldChangeTime, PMKtoLogiTask.get(key).tadHeur);// 第1个产品第1个机器第1个批次的物料的开始时间
        // 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            taf[prodSeq[p]] = timeCur;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
        maxTaf =Arrays.stream( taf).max().getAsDouble();
//		minTard=0;
        // 计算maxTard
//		double sum_processTime=0;
//		double max_processTime=-1;
//		for (LogiTask logiTask:logiTaskSet) {
//			sum_processTime=sum_processTime+logiTask.getProcessTime();
//			if (logiTask.getProcessTime()>max_processTime) {
//				max_processTime=logiTask.getProcessTime();
//			}
//		}
//		maxTard=sum_processTime/para.numVeh+(para.numVeh-1)*max_processTime/para.numVeh;

        // 计算maxTaf
//		prodObj=solveProdProblem();//子问题求解
//		logiObj=solveLogiProblem();//子问题求解
//		relaxObj = prodObj + logiObj;
//		String key;
//		double[] maxTaf=new double[para.numProd];
//		key = prodSeq[0] + "_" + 0 + "_" + 0;
//		timeCur = Math.max( PMKtoLogiTask.get(key).getTss(), para.moldChangeTime);// 第1个产品第1个机器第1个批次的物料的开始时间（生产问题求得），是整个系统的最早开始时间
//		for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//			for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//				if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//					timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)+para.getDistance(para.machSet.get(m-1), para.machSet.get(m))/para.lineSpeed;//+(para.jobshopLength/para.numMach)/para.lineSpeed;
//				}
//				for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//					key = prodSeq[p] + "_" + m + "_" + k;
//					// 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//					timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).getTad()); // 真正的开始时间是物料到达和最早开工时间取大；
//					timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//				}
//			}
//			maxTaf[p] = timeCur;
//			timeCur=timeCur-(para.numMach-1)*para.cycTime[prodSeq[p]]+para.moldChangeTime-para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach-1))/para.lineSpeed;//-para.jobshopLength/para.lineSpeed;
//		}
//		tafUB[para.numProd-1]=maxTaf[para.numProd-1];
//		for (int p = 0; p < para.numProd; p++) {
//			tafUB[p]=maxTaf[p];
//		}
        //////////////////
        // 计算makespan
//		for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//			makespan=makespan+para.prodSet.get(p).getCycleTime()*(para.numMach-1+para.demands[p]);
//		}
        tardLB = new double[para.numProd];
        double[] productionTimeWithoutLogiConstriants = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            productionTimeWithoutLogiConstriants[p] = para.prodSet.get(prodSeq[p]).getCycleTime()
                    * (para.numMach - 1 + para.demands[prodSeq[p]])
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        tafLB = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            tafLB[p] = Float.MAX_VALUE;
        }
        getTafLB();
//        tardLB[prodSeq[0]] = Math.max(tafLB[prodSeq[0]] - productionTimeWithoutLogiConstriants[prodSeq[0]], 0);
//        double benchTime = 0;
//        for (int p = 1; p < para.numProd; p++) {
//            benchTime = tardLB[prodSeq[p - 1]] + productionTimeWithoutLogiConstriants[prodSeq[p - 1]];
//            benchTime += productionTimeWithoutLogiConstriants[prodSeq[p]];
//            tardLB[prodSeq[p]] = Math.max(tafLB[prodSeq[p]] - benchTime, 0);
//        }


    }

    public static void main(String[] args) throws IOException, ClassNotFoundException, IloException {

        Instance para = new Instance(2, 2, 2, false, 1, 0.5, 3,
                "data\\demand2.txt");
        para.Initialization(1);
        SubgradMethod sgd = new SubgradMethod(para, 1, false);
        sgd.initialization();
//        int[] b={1,0};
//        double a = sgd.LocalSearch(b);
//        ArrayList<Integer> a=new ArrayList<Integer>();
//        a.add(1);
//        a.add(2);
//        ArrayList<Integer> b=sgd.listCopy(a);
//        Collections.swap(a,0,1);
//        System.out.println(a.get(0));
//        System.out.println(b.get(0));
    }

    public void getTafLB() {

        double[] minBatchPerMach = new double[para.numMach];
        ArrayList<Integer[]> batchOrderAse = new ArrayList<Integer[]>();
//		int[][] K={{2,3},{1,1}};
        for (int m = 0; m < para.numMach; m++) {
            Integer[] batchOrder_at_m = new Integer[para.numProd];
            for (int p = 0; p < para.numProd; p++) {
                batchOrder_at_m[p] = para.K[p][m];
            }
            Arrays.sort(batchOrder_at_m);
            batchOrderAse.add(batchOrder_at_m);
        }

        ArrayList<Integer[]> accumBatchOrder = new ArrayList<Integer[]>();
        for (int m = 0; m < para.numMach; m++) {
            Integer[] accu = new Integer[para.numProd];
            accu[0] = batchOrderAse.get(m)[0];
            for (int p = 1; p < para.numProd; p++) {
                accu[p] = accu[p - 1] + batchOrderAse.get(m)[p];
            }
            accumBatchOrder.add(accu);
        }

        double[] logiLB = new double[para.numProd];
        for (int num_p = 0; num_p < para.numProd; num_p++) {
            logiLB[num_p] = getLB_ParallelMachine(para, accumBatchOrder, num_p);
        }


//		for (int m = 0; m < para.numMach; m++) {
//			minBatchPerMach[m]=Integer.MAX_VALUE;
//		}
//		for (int p = 0; p < para.numProd; p++) {
//			for (int m = 0; m < para.numMach; m++) {
//				if (para.K[p][m]<minBatchPerMach[m]){
//					minBatchPerMach[m]=para.K[p][m];
//				}
//			}
//		}
//		System.out.println(minBatchPerMach);


//		double[] tss_final_of_p=new double[para.numProd];
        // 先生成生产顺序
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p,
                    1 / ((para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                            + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
                            / para.lineSpeed));
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());
        int[] prodSeq = new int[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
        }
        // 后物流顺序
        // 先得到tss,并记录对应的pmk
//		ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
        double[] tafWithoutLogiCons = new double[para.numProd];
        double[][][] tssWithoutLogiCons = new double[para.numProd][para.numMach][para.maxK];
        String key;
        double timeCur = 0;// earliestReadyTime;// 时间游标
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur
                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
                            / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//					TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
//					tssWithoutLogiConsList.add(tssWithoutLogiCons);
//								PMKtoLogiTask.get(key).setTss(timeCur);
                    tssWithoutLogiCons[prodSeq[p]][m][k] = timeCur;
                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        tafWithoutLogiCons[prodSeq[p]] = timeCur;
                    }
                }
            }
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
                    / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
//		Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
//		String[] pmk = new String[para.logiTaskSet.size()];
//		for (int s = 0; s < para.logiTaskSet.size(); s++) {
//			pmk[s] = tssWithoutLogiConsList.get(s).pmk;
//		}
//        double[] tafLB = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {
//            tafLB[p] = Float.MAX_VALUE;
//        }
//        double[] accTaf=new  double[para.numProd];
//        accTaf[0]=tafWithoutLogiCons[prodSeq[p]]
//        for (int p = 0; p < para.numProd; p++) {
//            accTaf[p]
//        }
        double[] minDelta= new double[para.numProd];
        int[] lastMach= new int[para.numProd];
        double minD=Float.MAX_VALUE;
        int lastM=0;
        for (int p = 0; p < para.numProd; p++) {
            minDelta[p]=Float.MAX_VALUE;
        }
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    double delta = tafWithoutLogiCons[prodSeq[p]] - tssWithoutLogiCons[prodSeq[p]][m][k];
                    if (delta<minDelta[prodSeq[p]]){
                        minDelta[prodSeq[p]]=delta;
                        lastMach[prodSeq[p]]=m;
                    }
                    if (delta<minD){
                        minD=delta;
                    }

//                    double tad = logiLB[p] - para.getDistance(para.depot, para.machSet.get(m)) / para.vehSpeed;
//                    double temp = tad + delta;
//                    if (temp < tafLB[p]) {
//                        tafLB[p] = temp;
//                    }
                }
            }
        }
        for (int p = 0; p < para.numProd; p++) {
            if (lastMach[p]>lastM){
                lastM=lastMach[p];
            }
        }

        double a=0;
//        double minDelta=0;//new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
////            minDelta[p]=Float.MAX_VALUE;
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
//                    double delta[p][][] = tafWithoutLogiCons[p] - tssWithoutLogiCons[p][m][k];
//                    if (delta<minDelta) {
//                        minDelta = delta;
//                    }
//                }
//            }
//        }
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            double delta = tafWithoutLogiCons[prodSeq[p]] - tssWithoutLogiCons[prodSeq[p]][m][k];
            double tad = logiLB[p] - para.getDistance(para.depot, para.machSet.get(lastM)) / para.vehSpeed;
            tafLB[p] = tad + minD;
        }
    }

    public double getLB_ParallelMachine(Instance para, ArrayList<Integer[]> accumBatchOrder, int num_p) {
        double LB = Float.MIN_VALUE;
        int[] num_task = new int[para.numMach];

        double[] processTime = new double[para.numMach];
        double totalProcTime = 0;
        for (int m = 0; m < para.numMach; m++) {
            num_task[m] = accumBatchOrder.get(m)[num_p];
            processTime[m] = (para.getDistance(para.depot, para.machSet.get(m)) + para.getDistance(para.machSet.get(m), para.depot)) / para.vehSpeed +
                    (para.getDistance(para.depot, para.machSet.get(m)) + para.getDistance(para.machSet.get(m), para.depot)) * para.consumingRate * para.chargingRate;
            totalProcTime += num_task[m] * processTime[m];
        }
        ArrayList<Double> allProcessTimes = new ArrayList<Double>();
        for (int m = 0; m < para.numMach; m++) {
            for (int i = 0; i < num_task[m]; i++) {
                allProcessTimes.add(processTime[m]);
            }
        }
        allProcessTimes.sort(Comparator.naturalOrder());
        double[] candidate_LB = new double[3];
        candidate_LB[0] = allProcessTimes.get(allProcessTimes.size() - 1);
//		if (allProcessTimes.size()-1-para.numMach>=0) {
        candidate_LB[1] = 0;//allProcessTimes.get(allProcessTimes.size() - 1 - para.numMach) + allProcessTimes.get(allProcessTimes.size() - 1 - para.numMach + 1);
//		}
        candidate_LB[2] = totalProcTime / para.numVeh;
        for (int i = 0; i < 3; i++) {
            if (candidate_LB[i] > LB) {
                LB = candidate_LB[i];
            }
        }
        for (int k = 1; k <= para.numVeh; k++) {
//			int alphaUB=(allProcessTimes.size()-k)/para.numMach;
            for (int l = 1; l <= allProcessTimes.size(); l++) {
//				int l=alpha*para.numMach+k;
                int lambda = k * (l / para.numVeh) + Math.min(k, l - (l / para.numVeh) * para.numVeh);
                ArrayList<Double> jobs = new ArrayList<Double>();
                for (int i = 0; i < lambda; i++) {
                    jobs.add(allProcessTimes.get(i));
                }
                double[] temp = new double[3];
                temp[0] = jobs.get(jobs.size() - 1);
                temp[1] = 0;
                double total = 0;
                for (int i = 0; i < jobs.size(); i++) {
                    total += jobs.get(i);
                }
                temp[2] = total / k;
                for (int i = 0; i < 3; i++) {
                    if (temp[i] > LB) {
                        LB = temp[i];
                    }
                }
            }
        }
        return LB;
    }

//    public void run3() throws ClassNotFoundException, IOException, IloException {
//        boolean flag = true;
//        bestLambdaForLowerBound = new double[logiTaskSet.size()];
//        double[] bestTssForLowerBound = new double[logiTaskSet.size()];
//        double[] bestTadForLowerBound = new double[logiTaskSet.size()];
//        // 为每一个logiTask设置其lambda（初始化拉格朗日乘子）
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
//                    PMKtoLogiTask.get(p + "_" + m + "_" + k).setLambda(Math.random()*0.);// (double)Math.random()*0.1
//                }
//            }
//        }
//        for (LogiTask logiTask : logiTaskSet) {
//            bestLambdaForLowerBound[logiTask.getId()] = 0.0;
//        }
//        iter = 0;
//        while (flag) {// (flag && iter<1000) {
//            iter++;
////			System.out.println("========iteration " + iter + ", incumbent bestRelaxObj="+bestRelaxObj+", bestObjUpBound="+bestObjUpBound +", bestIter="+bestIter+ " ============");
////			给定当前的lambda，求解lagrangian松弛问题 （最关键）
//            prodObj = solveProdProblemNoLift();// 子问题求解
//            logiObj = solveLogiProblem();// 子问题求解
//            LBNoLift = prodObj + logiObj;
//            LBLift =Math.max(LBLift,solveProdProblem()+logiObj);
////            if(LBLift>2310.0){
////                solveProdProblem();
////            }
//            // 给定prodSeq得到两个算法下的可行目标值即上界
//            double objUpBound1 = LagrangHeuristic();
//            if (objUpBound1 < bestObjUpBound) {
//                bestObjUpBound = objUpBound1;
//            }
//            UBNoImprove=Math.min(objUpBound1,UBNoImprove);
//
//            double temp = LocalSearch(prodSeq,logiSeq);
//            if (temp < bestObjUpBound) {
//                bestObjUpBound = temp;
//                for (int ii = 0; ii < prodSeq.length; ii++) {
//                    bestProdSeq[ii] = prodSeq[ii];
//                }
//
//                bestLogiSeq.clear();
//                for (int ii = 0; ii < logiSeq.size(); ii++) {
//                    ArrayList<Integer> alist = new ArrayList<Integer>();
//                    for (int jj = 0; jj < logiSeq.get(ii).size(); jj++) {
//                        alist.add(logiSeq.get(ii).get(jj));
//                    }
//                    bestLogiSeq.add(alist);
//                }
//            }
//            if (LBNoLift > bestRelaxObj) {// 本次迭代松弛问题的目标值有改进(提升)
//                improved = true;
//                bestRelaxObj = LBNoLift;
//                for (LogiTask logiTask : logiTaskSet) {
//                    bestLambdaForLowerBound[logiTask.getId()] = logiTask.getLambda();
//                    bestTssForLowerBound[logiTask.getId()] = logiTask.getTss();
//                    bestTadForLowerBound[logiTask.getId()] = logiTask.getTad();
//                }
//                iter_no_improve = 0;
//            } else {
//                improved = false;
//                iter_no_improve++;
//            }
////          获得梯度
//            updateSubgradient();// update lambda with dualSolAndObj.solution
////          更新步长
//            updateStepsize();
//            // 输出信息
//            if (show) {
//                InfoOutput();
//            }
////			更新乘子
//            updateLambda();
//            if (iter_no_improve > max_iter_no_improve) {
//                flag = false;
//            }
//            if (iter > 1 && bestObjUpBound <= LBLift) {
////                LBLift = bestObjUpBound;
//                flag = false;
//            }
//            if (iter == 1 && LBLift - bestObjUpBound > 0.005) {
//                System.out.println("wrong");
//            }
////            System.out.println("iter="+iter+", LB="+bestRelaxObj+", UB="+bestObjUpBound);
//        }
////        System.out.println(tafLB[para.numProd-1]<para.prodSet.get(prodSeq[para.numProd-1]).getTaf());
////        System.out.println(para.prodSet.get(prodSeq[para.numProd-1]).getTaf());
////        System.out.println(" ");
//    }

    public double ini_heur(int[] prodSeq,  ArrayList<ArrayList<Integer>> ini_logiSeq) {

        double bestObjVal=Float.MAX_VALUE;

        for (int ii = 0; ii < 2000; ii++) {

            if (Math.random() < 0.5) {
                Random rnd = new Random();
                int p_change = rnd.nextInt(para.numProd - 1);
                int temp = prodSeq[p_change + 1];
                prodSeq[p_change + 1] = prodSeq[p_change];
                prodSeq[p_change] = temp;
            }
            ArrayList<ArrayList<Integer>> cumb_logiSeq = new ArrayList<ArrayList<Integer>>();
            for (int r = 0; r < para.numVeh; r++) {
                cumb_logiSeq.add(new ArrayList<Integer>());
            }

            ArrayList<Integer> first = new ArrayList<>();
            ArrayList<Integer> last = new ArrayList<>();
            for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
                int pindex = prodSeq[p];
                ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
                ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
                for (int index = 0; index < arrayofPMKI.size(); index++) {
                    arrayOfI.add(arrayofPMKI.get(index)[3]);
                }
                first.add(arrayOfI.get(0));
                last.add(arrayOfI.get(arrayOfI.size() - 1));
            }

            ArrayList<Integer> prodPriority = new ArrayList<Integer>();
            for (int p = 0; p < prodSeq.length; p++) {
                prodPriority.add(prodSeq[p]);
            }
            int[] priorityLogiTask = new int[logiTaskSet.size()];
            for (int i = 0; i < ini_logiSeq.size(); i++) {
                for (int j = 0; j < ini_logiSeq.get(i).size(); j++) {
                    int logiTaskId = ini_logiSeq.get(i).get(j);
                    int prodId = para.ItoLogiTask.get(logiTaskId).getProd().getId();
                    priorityLogiTask[logiTaskId] = prodPriority.indexOf(prodId);
                }
            }
            for (int i = 0; i < ini_logiSeq.size(); i++) {
                ArrayList<LogiTaskPrior> priorSet = new ArrayList<LogiTaskPrior>();
                for (int j = 0; j < ini_logiSeq.get(i).size(); j++) {
                    int logiTaskId = ini_logiSeq.get(i).get(j);
                    priorSet.add(new LogiTaskPrior(logiTaskId, priorityLogiTask[logiTaskId]));
//                System.out.print(priorityLogiTask[logiTaskId]+" ");
                }
//            System.out.println(priorSet.size());
                Collections.sort(priorSet, new LogiTaskSortByPrior());

                for (int j = 0; j < priorSet.size(); j++) {
                    cumb_logiSeq.get(i).add(priorSet.get(j).id);
                }

            }

            bestObjVal = Math.min(getObj(prodSeq, cumb_logiSeq), getObj(prodSeq, ini_logiSeq));

            ArrayList<ArrayList<Integer>> logiSeq = listCopy2(cumb_logiSeq);

            int iter = 0;
            Random rnd = new Random();
            while (iter < 20) {
                int ov;
                int op;
                while (true) {
                    int a = rnd.nextInt(para.numVeh);
                    if (logiSeq.get(a).size() > 1) {
                        ov = a;
                        op = rnd.nextInt(logiSeq.get(a).size());
                        break;
                    }
                }
                int o_LogiTaskId = logiSeq.get(ov).get(op);
                int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
                int o_Priority = prodPriority.indexOf(o_ProdId);

                HashSet<int[]> feasInnerInsertPos = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                HashSet<int[]> feasAfterInsertPos = new HashSet<>();//插入在指定位置之后
                HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
                HashSet<int[]> feasSwapPos = new HashSet<>();//相同优先级的交换

                for (int i = 0; i < logiSeq.size(); i++) {
                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
                        int[] fea = {i, j};
                        int n_LogiTaskId = logiSeq.get(i).get(j);
                        if (o_LogiTaskId != n_LogiTaskId) {
                            int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                            int n_Priority = prodPriority.indexOf(n_ProdId);
                            if (n_Priority == o_Priority) {
                                feasSwapPos.add(fea);
                                feasInnerInsertPos.add(fea);
                            }
                            if (j < logiSeq.get(i).size() - 1) {
                                int n2_LogiTaskId = logiSeq.get(i).get(j + 1);
                                int n2_ProdId = para.ItoLogiTask.get(n2_LogiTaskId).getProd().getId();
                                int n2_Priority = prodPriority.indexOf(n2_ProdId);
                                if (n_Priority <= o_Priority && n2_Priority > o_Priority) {
                                    feasAfterInsertPos.add(fea);
                                }
                            } else {
                                if (n_Priority <= o_Priority) {
                                    feasEndInsertPos.add(fea);
                                }
                            }
                        }
                    }
                }
                Random rndOperator = new Random();
                boolean improved = false;
                while (!improved && (!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty() || !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty())) {
                    switch (rndOperator.nextInt(4)) {
                        case 0:
                            if (!feasInnerInsertPos.isEmpty()) {
                                Iterator<int[]> iterator = feasInnerInsertPos.iterator();
                                int[] fea = iterator.next();
                                feasInnerInsertPos.remove(fea);
                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                                innerInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
                                double testObjVal = getObj(prodSeq, neigh_logiSeq);
                                if (testObjVal <= bestObjVal) {
                                    bestObjVal = testObjVal;
                                    logiSeq = neigh_logiSeq;
//                                if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                    int a=0;
//                                }
                                    improved = true;
                                }
                            }
                            break;
                        case 1:
                            if (!feasEndInsertPos.isEmpty()) {
                                Iterator<int[]> iterator = feasEndInsertPos.iterator();
                                int[] fea = iterator.next();
                                feasEndInsertPos.remove(fea);
                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                                endInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                                double testObjVal = getObj(prodSeq, neigh_logiSeq);
                                if (testObjVal <= bestObjVal) {
                                    bestObjVal = testObjVal;
                                    logiSeq = neigh_logiSeq;
                                    improved = true;
                                }
                            }
                            break;
                        case 2:
                            if (!feasAfterInsertPos.isEmpty()) {
                                Iterator<int[]> iterator = feasAfterInsertPos.iterator();
                                int[] fea = iterator.next();
                                feasAfterInsertPos.remove(fea);
                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                                afterInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                                double testObjVal = getObj(prodSeq, neigh_logiSeq);
                                if (testObjVal <= bestObjVal) {
                                    bestObjVal = testObjVal;
                                    logiSeq = neigh_logiSeq;
                                    improved = true;
                                }
                            }
                            break;
                        case 3:
                            if (!feasSwapPos.isEmpty()) {
                                Iterator<int[]> iterator = feasSwapPos.iterator();
                                int[] fea = iterator.next();
                                feasSwapPos.remove(fea);
                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                                swapOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                                double testObjVal = getObj(prodSeq, neigh_logiSeq);
                                if (testObjVal <= bestObjVal) {
                                    bestObjVal = testObjVal;
                                    logiSeq = neigh_logiSeq;
                                    improved = true;
                                }
                            }
                            break;
                    }
                }
                iter++;
            }

            int[][] posi = new int[para.numProd][4];
            for (int i = 0; i < logiSeq.size(); i++) {
                for (int j = 0; j < logiSeq.get(i).size(); j++) {
                    int logiTaskId = logiSeq.get(i).get(j);
                    int prodId = para.ItoLogiTask.get(logiTaskId).getProd().getId();
                    if (first.contains(logiTaskId)) {
                        posi[prodId][0] = i;
                        posi[prodId][1] = j;
                    }
                    if (last.contains(logiTaskId)) {
                        posi[prodId][2] = i;
                        posi[prodId][3] = j;
                    }
                }
            }
            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
            bestObjVal = Math.min(bestObjVal, LocalImprove(prodSeq, neigh_logiSeq, bestObjVal, posi));
        }

        return bestObjVal;
    }

    public int[] getIniProdSeq(){
        int[] prodSeq=new int[para.numProd];

        double[] sum_lambda_wrt_prod = new double[para.numProd];
        double[] denominator_wrt_prod = new double[para.numProd];

        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
            int corresProdId = logiTask.getProd().getId();
            sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + logiTask.getLambda();
        }
        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p,
                    (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列

        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
        }
        return prodSeq;
    }
    public ArrayList<ArrayList<Integer>> getIniLogiSeq(int[] prodSeq){

        //////////////////////////
        //////////////////////////

        double[] timeCur = new double[para.numVeh];
        ArrayList<ArrayList<Integer>> logiSeq3 = new ArrayList<ArrayList<Integer>>();
        for (int r = 0; r < para.numVeh; r++) {
//			preLogiTaskID[r]=-1;
            logiSeq3.add(new ArrayList<Integer>());
        }
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            ArrayList<Integer> idSet = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }

            while (!arrayOfI.isEmpty()) {
                Random random = new Random();
                int index = random.nextInt(arrayOfI.size());
                idSet.add(arrayOfI.get(index));
                arrayOfI.remove(index);
            }

            for (int index = 0; index < idSet.size(); index++) {
//					int machId = arrayofPMKI.get(index)[1];// 对应mach的ID
//					int batchId = arrayofPMKI.get(index)[2];// 对应batch的ID
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
//						Random random2 = new Random();
//						int index2=random2.nextInt(para.numVeh);
//						logiSeq.get(index2).add(logiTaskId);
                int min_index = 0;
                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                double min_val = Float.MAX_VALUE;
                for (int r = 0; r < para.numVeh; r++) {
                    if (timeCur[r] + requiredDeliTime < min_val) {
                        min_val = timeCur[r] + requiredDeliTime;
                        min_index = r;
                    }
                }
                logiSeq3.get(min_index).add(logiTaskId);
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }

        }
        /////////////////////////////
        /////////////////////////////

        return logiSeq3;
    }

    public double[] arrayCopy(double[] array){
        double[] a2=new double[array.length];
        for (int i = 0; i < array.length; i++) {
            a2[i]=array[i];
        }
        return a2;
    }

    public void run_surrogate(double z_FCFS, int[] f_prodSeq, ArrayList<ArrayList<Integer>> f_logiSeq) throws ClassNotFoundException, IOException, IloException {
        boolean flag = true;
        double iniObj=ini_heur(f_prodSeq,f_logiSeq);

        double[] tad_sur_cur=new double[para.logiTaskSet.size()];
        double[] tss_sur_cur=new double[para.logiTaskSet.size()];
        double[] lambda_sur_cur=new double[para.logiTaskSet.size()];
        double[] taf_sur_cur=new double[para.numProd];
        double[] tad_sur_last=new double[para.logiTaskSet.size()];
        double[] tss_sur_last=new double[para.logiTaskSet.size()];
        double[] taf_sur_last=new double[para.prodSet.size()];
        double[] lambda_sur_last=new double[para.logiTaskSet.size()];

        // 为每一个logiTask设置其lambda（初始化拉格朗日乘子）
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
                    lambda_sur_cur[PMKtoLogiTask.get(p + "_" + m + "_" + k).getId()]=0.0;
                }
            }
        }

        int[] prodSeq_ini=getIniProdSeq();
//        printArray("prodSeq",prodSeq_ini);

        ArrayList<ArrayList<Integer>> logiSeq_ini=getIniLogiSeq(prodSeq_ini);
//        printlogiSeq(logiSeq_ini);

        iniObj=Math.min(ini_heur(prodSeq_ini,logiSeq_ini),iniObj);
//        printNumber("qhat",iniObj);

        bestObjUpBound = Math.min(z_FCFS,iniObj);//LagrangHeuristic();
//        printNumber("bestUB= min(qhat,z_FCFS)",bestObjUpBound);

        int iter=1;
        double M=10;
        double r=0.2;
        double surDualEstimate=iniObj;
        double p,alpha;

        int[] prodSeq_sur_cur=getSurProdSeq(lambda_sur_cur);
        ArrayList<ArrayList<Integer>> logiSeq_sur_cur=getSurLogiSeq_ini(lambda_sur_cur);

        tad_sur_cur=getSurTad(logiSeq_sur_cur);
        taf_sur_cur=getSurTaf(prodSeq_sur_cur,lambda_sur_cur);
        tss_sur_cur=getSurTss(taf_sur_cur,prodSeq_sur_cur);

        double surDual = getSurProdObj(lambda_sur_cur,tss_sur_cur,taf_sur_cur) + getSurLogiObj(lambda_sur_cur,tad_sur_cur);
//        printNumber("surDual",surDual);

        bestRelaxObj=Math.max(surDual,bestRelaxObj);
//        printNumber("bestLB",bestRelaxObj);

        double[] surSubGrad=getSurGrad(tad_sur_cur,tss_sur_cur);
        double stepsize=(surDualEstimate-surDual)/Math.pow(getSurNorm(surSubGrad),2);
//        System.out.println("c0:"+stepsize);
        tad_sur_last=arrayCopy(tad_sur_cur);
//                printArray("tad_sur_last",tad_sur_last);
        tss_sur_last=arrayCopy(tss_sur_cur);
//                printArray("tss_sur_last",tss_sur_last);
        taf_sur_last=arrayCopy(taf_sur_cur);
//                printArray("taf_sur_last",taf_sur_last);

        while (flag) {// (flag && iter<1000) {
//            System.out.println("============================");
//            System.out.print("iteration "+iter+"  ");
            p = 1-1/Math.pow(iter,r);
            alpha = 1-1/( M*Math.pow(iter,p) );
//            printNumber("alpha",alpha);
//
//            if(iter==0){
//
//            }

            double[] surGrad_last=getSurGrad(tad_sur_last,tss_sur_last);
//            printArray("grad_sur_last",surGrad_last);

            double[] surGrad_cur=getSurGrad(tad_sur_cur,tss_sur_cur);
//            printArray("grad_sur_current",surGrad_cur);

            double norm_surGrad_last=getSurNorm(surGrad_last);
            double norm_surGrad_cur=getSurNorm(surGrad_cur);
//            System.out.println("ratio="+norm_surGrad_last/norm_surGrad_cur);

            stepsize=alpha*stepsize*norm_surGrad_last/norm_surGrad_cur;
//            System.out.println("c"+iter+"="+stepsize);
//            lambda_sur_last=arrayCopy(lambda_sur);
            lambda_sur_last=arrayCopy(lambda_sur_cur);
//            printArrayNoFormat("lambda_sur_last",lambda_sur_last);

            lambda_sur_cur=updateSurLambda(surGrad_cur,lambda_sur_cur,stepsize);
//            printArrayNoFormat("lambda_sur_current",lambda_sur_cur);



            tad_sur_last=arrayCopy(tad_sur_cur);
//            printArray("tad_sur_last",tad_sur_last);
            tss_sur_last=arrayCopy(tss_sur_cur);
//            printArray("tss_sur_last",tss_sur_last);
            taf_sur_last=arrayCopy(taf_sur_cur);
//            printArray("taf_sur_last",taf_sur_last);

            double surDualTarget=getSurDual(lambda_sur_cur,tad_sur_last,tss_sur_last,taf_sur_last);
//            printNumber("\033[31msurDualTarget\033[0m",surDualTarget);

            double prodObj_using_last_x=getSurProdObj(lambda_sur_cur,tss_sur_last,taf_sur_last);
//            printNumber("\033[31mprodObj using last x\033[0m",prodObj_using_last_x);

            prodSeq_sur_cur=getSurProdSeq(lambda_sur_cur);
            taf_sur_cur=getSurTaf(prodSeq_sur_cur,lambda_sur_cur);
//            printArray("taf_sur",taf_sur_cur);

            tss_sur_cur=getSurTss(taf_sur_cur,prodSeq_sur_cur);
//            printArray("tss_sur",tss_sur_cur);

            double prodObj_cur=getSurProdObj(lambda_sur_cur,tss_sur_cur,taf_sur_cur);
//            printNumber("\033[31mprodObj using new x\033[0m",prodObj_cur);

//            if(prodObj_using_last_x<prodObj_cur){
//                System.out.println("\033[31mwrong: prodObj_last <= prodObj_new\033[0m"+ "prodLast="+prodObj_using_last_x+",  prodCur="+prodObj_cur);
////                double[] a = getSurTaf(prodSeq_sur_cur,lambda_sur_cur);
//            }
//            printNumber("new prodObj",prodObj_cur);

//            printNumber("\033[31mnew logiObj should be less than \033[0m", surDualTarget-prodObj_cur);

            logiSeq_sur_cur=getSurLogiSeq_Heur(lambda_sur_cur,surDualTarget-prodObj_cur,logiSeq_sur_cur);

            tad_sur_cur=getSurTad(logiSeq_sur_cur);
            double logiObj_cur=getLogiObjSur(lambda_sur_cur,tad_sur_cur);
//            printNumber("\033[31mnew logiObj\033[0m: ",logiObj_cur);

            surDual=prodObj_cur+logiObj_cur;
//            printNumber("\033[31mnew surDual\033[0m",surDual);

//            printArray("tad_sur",tad_sur_cur);

//            printNumber("new logiObj",logiObj_cur);

//            if(logiObj_cur>surDualTarget-prodObj_cur){
//                System.out.println("\033[31mwrong: logiObj_new + prodObj_new> targetSurDual\033[0m"+ "logiTarget="+(surDualTarget-prodObj_cur)+",  logiCur="+logiObj_cur);
////                double a=getSurDual(lambda_sur_cur,tad_sur_last,tss_sur_last,taf_sur_last);
////                System.out.println();
//            }

//            printArray("prodSeq",prodSeq_sur_cur);
//            printlogiSeq(logiSeq_sur_cur);

            bestRelaxObj=surDual;//Math.max(surDual,bestRelaxObj);
            double liftedProdObj=liftTaf(tafLB,prodSeq_sur_cur,lambda_sur_cur);
//            LBLift=Math.max(bestRelaxObj,LBLift);
            LBLift=Math.max(liftedProdObj+logiObj_cur,LBLift);


            double temp = LocalSearch(prodSeq_sur_cur,logiSeq_sur_cur);
            UBNoImprove=bestObjUpBound;
            if (temp < bestObjUpBound) {
                bestObjUpBound = temp;
            }
            

            double[] diff_lambda_sur=new double[para.logiTaskSet.size()];
            for (int i = 0; i < diff_lambda_sur.length; i++) {
                diff_lambda_sur[i]=lambda_sur_cur[i]-lambda_sur_last[i];
            }
            double diff=getSurNorm(diff_lambda_sur);
//            printNumber("diff_lambda_norm",diff);
            if(diff<0.1){
                flag=false;
            }
            iter++;
        }

    }

    private double liftTaf( double[] tafLB, int[] prodSeq,double[] lambda) {
        double liftedProdObj=0;

        String key;
        double[] taf=new double[para.prodSet.size()];
        double timeCur =minTard;//Math.max(minTard,tardLB[prodSeq[0]]);// minTard;//para.moldChangeTime;//earliestReadyTime;// 时间游标
        double[] tempTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        timeCur=Math.max( timeCur, tafLB[p] );
                        tempTaf[prodSeq[p]] = timeCur;
                    }
                }
            }
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }


        double[] latestTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            latestTaf[p] = tempTaf[p] + maxTaf - Arrays.stream(tempTaf).max().getAsDouble();//tempTaf[prodSeq[para.numProd - 1]];
        }

        double[] sum_lambda_wrt_prod = new double[para.numProd];
        double[] denominator_wrt_prod = new double[para.numProd];
        for (int i = 0; i < sum_lambda_wrt_prod.length; i++) {
            sum_lambda_wrt_prod[i]=0;
            denominator_wrt_prod[i]=0;
        }
        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
            int corresProdId = logiTask.getProd().getId();
            int corresLogiId = logiTask.getId();
            sum_lambda_wrt_prod[corresProdId] +=  lambda[corresLogiId];
        }
        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }

        double[] ratio=new double[para.numProd];

        for (int p = 0; p < para.numProd; p++) {
            ratio[p]=(1.0 - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p];
        }

        for (int p = 0; p < para.numProd; p++) {
            if (ratio[p] < 0)//只有总权重小于0才会令其Taf尽可能大
            {
                taf[p]=latestTaf[p];
            }else {
                taf[p]=tempTaf[p];
            }
        }

        // 算最终的目标函数
        liftedProdObj = 0;
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            liftedProdObj +=  taf[prodSeq[p]];
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    String key2 = prodSeq[p] + "_" + m + "_" + k;
                    double tss = taf[prodSeq[p]]
                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
                            / para.lineSpeed;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    liftedProdObj -=  tss * lambda[ PMKtoLogiTask.get(key2).getId()];
                }
            }
        }


        return liftedProdObj;
    }

    private double getLogiObjSur(double[] lambda, double[] tad) {
        double obj = 0 ;
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            obj += lambda[i] * tad[i] ;
        }
        return obj;
    }

    private double[] getSurTad(ArrayList<ArrayList<Integer>> logiSeq) {
        double[] tad=new double[para.logiTaskSet.size()];

        for (int i = 0; i < logiSeq.size(); i++) {
            double baseTime = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {
                int logiTaskId = logiSeq.get(i).get(j);
                tad[logiTaskId]=baseTime + para.ItoLogiTask.get(logiTaskId).getChargingTime() + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2 ;
                baseTime = baseTime + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }
        }

        return tad;
    }

    private ArrayList<ArrayList<Integer>> getSurLogiSeq_ini(double[] lambda) {

        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<>();
        for (int r = 0; r < para.numVeh; r++) {
            logiSeq.add(new ArrayList<Integer>());
        }


        ArrayList<LogiTask> logiTaskList=new ArrayList<>();
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            LogiTask logi=new LogiTask();
            logi.setId(i);
            logi.setLambda(lambda[i]);
            logi.setProcessingTime(para.ItoLogiTask.get(i).getProcessTime());
            logi.setTransTime(para.ItoLogiTask.get(i).getTransTime());
            logi.setChargingTime(para.ItoLogiTask.get(i).getChargingTime());
            logiTaskList.add(logi);
        }

        Collections.sort(logiTaskList,new LogiTaskSortByLambdaProcessTimeRatio());
//        for (int i = 0; i < para.numVeh; i++) {
//            System.out.println(checkLambdaProcessTimeRatio1(logiTaskList));
//        }

        double[] timeCur=new double[para.numVeh];
        for (int i = 0; i < para.numVeh; i++) {
            timeCur[i]=0;
        }
        int min_index=-1;
        for (int i=0;i<logiTaskList.size();i++) {
            LogiTask logi=logiTaskList.get(i);
            double requiredDeliTime = logi.getChargingTime()
                    + logi.getTransTime() / 2;
            double min_val = Float.MAX_VALUE;
            for (int r = 0; r < para.numVeh; r++) {
                if (timeCur[r] + requiredDeliTime < min_val) {
                    min_val = timeCur[r] + requiredDeliTime;
                    min_index = r;
                }
            }
            logiSeq.get(min_index).add(logi.getId());
            timeCur[min_index] = timeCur[min_index] + logi.getProcessTime();
        }

        return logiSeq;
    }

    private double getSurProdObj(double[] lambda, double[] tss, double[] taf) {
        double obj=0;
        for (int i = 0; i < para.numProd; i++) {
            obj+=taf[i];
        }
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            obj-=lambda[i]*tss[i];
        }
        return obj;
    }

    private double getSurLogiObj(double[] lambda, double[] tad) {
        double obj=0;
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            obj+=lambda[i]*tad[i];
        }
        return obj;
    }

    private double[] getSurTaf(int[] prodSeq, double[] lambda) {
        double[] taf=new double[para.numProd];

        String key;
        double timeCur = 0; // 时间游标
        double[] tempTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        tempTaf[prodSeq[p]] = timeCur;
                    }
                }
            }
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }

        double[] latestTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            latestTaf[p] = tempTaf[p] + maxTaf - Arrays.stream(tempTaf).max().getAsDouble();//tempTaf[prodSeq[para.numProd - 1]];
        }

        double[] sum_lambda_wrt_prod = new double[para.numProd];
        double[] denominator_wrt_prod = new double[para.numProd];
        for (int i = 0; i < sum_lambda_wrt_prod.length; i++) {
            sum_lambda_wrt_prod[i]=0;
            denominator_wrt_prod[i]=0;
        }
        
        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
            int corresProdId = logiTask.getProd().getId();
            int corresLogiId = logiTask.getId();
            sum_lambda_wrt_prod[corresProdId] +=  lambda[corresLogiId];
        }
        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        double[] ratio=new double[para.numProd];

        for (int p = 0; p < para.numProd; p++) {
            ratio[p]=(1.0 - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p];
        }

        for (int p = 0; p < para.numProd; p++) {
            if (ratio[p] < 0)//只有总权重小于0才会令其Taf尽可能大
            {
                taf[p]=latestTaf[p];
            }else {
                taf[p]=tempTaf[p];
            }
        }

        return taf;
    }

    private double[] getSurTss(double[] taf,int[] prodSeq) {
        double[] tss=new double[para.logiTaskSet.size()];
        // 基于taf更新tss
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    double tss_temp = taf[prodSeq[p]]
                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
                            / para.lineSpeed;
                    tss[PMKtoLogiTask.get(prodSeq[p] + "_" + m + "_" + k).getId()]=tss_temp;
                }
            }
        }
        return tss;
    }

    private int[] getSurProdSeq(double[] lambdaSur) {
        int[] prodSeq=new int[para.numProd];
        double[] sum_lambda_wrt_prod = new double[para.numProd];
        double[] denominator_wrt_prod = new double[para.numProd];
        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
            int corresProdId = logiTask.getProd().getId();
            sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + lambdaSur[corresProdId];
        }
        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p,
                    (1 - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列


        for (int p = 0; p < para.numProd; p++) {
//            prodSeq[prodList.get(p).getId()] = p;//prodList.get(p).getId();
            prodSeq[p] = prodList.get(p).getId() ;
        }

        return prodSeq;
    }

//    public void run_surrogate_back(double z_FCFS, int[] f_prodSeq, ArrayList<ArrayList<Integer>> f_logiSeq) throws ClassNotFoundException, IOException, IloException {
//        boolean flag = true;
//        double iniObj=ini_heur(f_prodSeq,f_logiSeq);
//        double[] tad_sur=new double[para.logiTaskSet.size()];
//        double[] tss_sur=new double[para.logiTaskSet.size()];
//        double[] lambda_sur=new double[para.logiTaskSet.size()];
//        double[] taf_sur=new double[para.numProd];
//        // 为每一个logiTask设置其lambda（初始化拉格朗日乘子）
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
//                    PMKtoLogiTask.get(p + "_" + m + "_" + k).setLambda(Math.random()*0.);// (double)Math.random()*0.1
//                    lambda_sur[PMKtoLogiTask.get(p + "_" + m + "_" + k).getId()]=PMKtoLogiTask.get(p + "_" + m + "_" + k).getLambda();
//                }
//            }
//        }
//
//
//        int[] prodSeq=getIniProdSeq();
//        printArray("prodSeq",prodSeq);
//
//        ArrayList<ArrayList<Integer>> logiSeq=getIniLogiSeq(prodSeq);
//        printlogiSeq(logiSeq);
//
//        iniObj=Math.min(ini_heur(prodSeq,logiSeq),iniObj);
//        printNumber("qhat",iniObj);
////        for (int i = 0; i < logiTaskSet.size(); i++) {
////            System.out.println("id="+logiTaskSet.get(i).getId()+" ,tad="+logiTaskSet.get(i).getTad()+", tss="+logiTaskSet.get(i).getTss()+" ,grad="+logiTaskSet.get(i).getSubgrad()+", lambda="+logiTaskSet.get(i).getLambda());
////        }
////        System.out.println();
//
//        bestObjUpBound = Math.min(z_FCFS,iniObj);//LagrangHeuristic();
//        printNumber("bestUB= min(qhat,z_FCFS)",bestObjUpBound);
//
//        iter=0;
//        double M=25;
//        double r=0.06;
//        double surDualEstimate=iniObj;
//        double p,alpha;
//
//        double[] tad_sur_last=new double[para.logiTaskSet.size()];
//        double[] tss_sur_last=new double[para.logiTaskSet.size()];
//        double[] taf_sur_last=new double[para.prodSet.size()];
//        double[] lambda_sur_last=new double[para.logiTaskSet.size()];
//
//        double surDual = solveProdProblemNoLift(lambda_sur) + solveLogiProblem_Ini(lambda_sur);//得到初始的tad_sur和tss_sur以及taf_sur
//        printNumber("surDual",surDual);
//
//        bestRelaxObj=Math.max(surDual,bestRelaxObj);
//        printNumber("bestLB",bestRelaxObj);
//
//        double denominator = 0;
//        for (LogiTask logiTask : logiTaskSet) {
//            denominator = denominator
//                    + ((logiTask.getTad() - logiTask.getTss()) * (logiTask.getTad() - logiTask.getTss()));
//        }
//        stepsize=(surDualEstimate-surDual)/denominator;
//        System.out.println("c0:"+stepsize);
//
//        while (flag) {// (flag && iter<1000) {
//            System.out.println("============================");
//            System.out.println("iteration "+iter);
//
//            p = 1-1/Math.pow(iter,r);
//            alpha = 1-1/( M*Math.pow(iter,p) );
//            printNumber("alpha",alpha);
//
//            if(iter==0){
//                tad_sur_last=arrayCopy(tad_sur);
//                printArray("tad_sur_last",tad_sur_last);
//                tss_sur_last=arrayCopy(tss_sur);
//                printArray("tss_sur_last",tss_sur_last);
//                taf_sur_last=arrayCopy(taf_sur);
//                printArray("taf_sur_last",taf_sur_last);
//            }
//
//            double[] surGrad_last=getSurGrad(tad_sur_last,tss_sur_last);
//            printArray("grad_sur_last",surGrad_last);
//
//            double[] surGrad=getSurGrad(tad_sur,tss_sur);
//            printArray("grad_sur_current",surGrad);
//
//            double norm_surGrad_last=getSurNorm(surGrad_last);
//            double norm_surGrad=getSurNorm(surGrad);
//            System.out.println("ratio="+norm_surGrad_last/norm_surGrad);
//
//            stepsize=alpha*stepsize*norm_surGrad_last/norm_surGrad;
//            System.out.println("c"+iter+"="+stepsize);
////            lambda_sur_last=arrayCopy(lambda_sur);
//            lambda_sur_last=arrayCopy(lambda_sur);
//            printArrayNoFormat("lambda_sur_last",lambda_sur_last);
//
//            lambda_sur=updateSurLambda(surGrad,lambda_sur,stepsize);
//            printArrayNoFormat("lambda_sur_current",lambda_sur);
//
//            double surDualTarget=getSurDual(lambda_sur,tad_sur_last,tss_sur_last,taf_sur_last);
//            printNumber("surDualTarget",surDualTarget);
//
//            tad_sur_last=arrayCopy(tad_sur);
//            printArray("tad_sur_last",tad_sur_last);
//            tss_sur_last=arrayCopy(tss_sur);
//            printArray("tss_sur_last",tss_sur_last);
//            taf_sur_last=arrayCopy(taf_sur);
//            printArray("taf_sur_last",taf_sur_last);
//
//            double prodObj_last=solveProdProblemNoLift(lambda_sur_last);
//            double prodObj=solveProdProblemNoLift(lambda_sur);
//            if(prodObj_last<=prodObj){
//                System.out.println("wrong: prodObj_last <= prodObj");
//                double a_last=solveProdProblemNoLift(lambda_sur_last);
//                double a=solveProdProblemNoLift(lambda_sur);
//                System.out.println();
//            }
//            printNumber("new prodObj",prodObj);
//            double logiObj=solveLogiProblem_Heur(lambda_sur,surDualTarget-prodObj);
//
//            surDual=prodObj+logiObj;
//            printNumber("new surDual",surDual);
//
//            printArray("tad_sur",tad_sur);
//            printArray("tss_sur",tss_sur);
//            printArray("taf_sur",taf_sur);
//            printNumber("new logiObj",logiObj);
//
//            printArray("prodSeq",prodSeq);
//            printlogiSeq(logiSeq);
//
//            bestRelaxObj=Math.max(surDual,bestRelaxObj);
//            LBLift=Math.max(bestRelaxObj,LBLift);
//
//            double temp = LocalSearch(this.prodSeq,this.logiSeq);
//            if (temp < bestObjUpBound) {
//                bestObjUpBound = temp;
//            }
//            UBNoImprove=bestObjUpBound;
//
//            double[] diff_lambda_sur=new double[para.logiTaskSet.size()];
//            for (int i = 0; i < diff_lambda_sur.length; i++) {
//                diff_lambda_sur[i]=lambda_sur[i]-lambda_sur_last[i];
//            }
//            double diff=getSurNorm(diff_lambda_sur);
//            printNumber("diff_lambda_norm",diff);
//            if(diff<0.00001){
//                flag=false;
//            }
//
////            updateSubgradient();//       获得梯度
//////            for (int i = 0; i < logiTaskSet.size(); i++) {
//////                System.out.println("id="+logiTaskSet.get(i).getId()+" ,tad="+logiTaskSet.get(i).getTad()+", tss="+logiTaskSet.get(i).getTss()+" ,grad="+logiTaskSet.get(i).getSubgrad()+", lambda="+logiTaskSet.get(i).getLambda());
//////            }
//////            System.out.println();
////
//////            System.out.println("before stepsize="+stepsize);
////            updateStepsize(surDual);//          更新步长
//////            System.out.println("after stepsize="+stepsize);
//////            System.out.println();
////
////            updateLambda();//			更新乘子
////            for (int i = 0; i < logiTaskSet.size(); i++) {
////                System.out.println("id="+logiTaskSet.get(i).getId()+" ,tad="+logiTaskSet.get(i).getTad()+", tss="+logiTaskSet.get(i).getTss()+" ,grad="+logiTaskSet.get(i).getSubgrad()+", lambda="+logiTaskSet.get(i).getLambda());
////            }
////
////            double targetSurDual=getSurDual( lambda_sur,tad_sur,tss_sur,taf_sur);//更新lambda后用现有的x带入得到
////
////            prodObj = solveProdProblemNoLift();// 子问题求解
////            logiObj = solveLogiProblem_Heur(lambda_sur,targetSurDual - prodObj );// 更新lambda后找一个比现有x更好的x
////            surDual=prodObj + logiObj;
////
////            LBNoLift = surDual;
////            LBLift =Math.max(LBLift,solveProdProblem()+logiObj);
////
////            double objUpBound1 = LagrangHeuristic();
////            if (objUpBound1 < bestObjUpBound) {
////                bestObjUpBound = objUpBound1;
////            }
////            UBNoImprove=Math.min( objUpBound1, UBNoImprove );
////
////            double temp = LocalSearch( prodSeq, logiSeq );
////            if (temp < bestObjUpBound) {
////                bestObjUpBound = temp;
////            }
////            if (LBNoLift > bestRelaxObj) {// 本次迭代松弛问题的目标值有改进(提升)
////                improved = true;
////                bestRelaxObj = LBNoLift;
////                iter_no_improve = 0;
////            } else {
////                improved = false;
////                iter_no_improve++;
////            }
////
////            if (iter_no_improve > max_iter_no_improve) {
////                flag = false;
////            }
////            if (bestObjUpBound <= LBLift) {
////                flag = false;
////            }
//            iter++;
//        }
//
//    }

    private void printNumber(String name, double number) {
        System.out.print(name);
        System.out.printf("%14.2f",number);
        System.out.println( );
    }

    private void printArray(String name, int[] array) {
        System.out.print(name+": ");
        for (int i = 0; i < array.length; i++) {
            System.out.print(array[i]+" ");
        }
        System.out.println( );
    }

    private void printArray(String name, double[] array) {
        System.out.print(name+": ");
        for (int i = 0; i < array.length; i++) {
            System.out.printf("%14.2f",array[i]);
            System.out.print(", ");
        }
        System.out.println( );
    }

    private void printArrayNoFormat(String name, double[] array) {
        System.out.print(name+": ");
        for (int i = 0; i < array.length; i++) {
            System.out.print(array[i]);
            System.out.print(", ");
        }
        System.out.println( );
    }

    private void printlogiSeq(ArrayList<ArrayList<Integer>> logiSeq) {
        System.out.println("logiSeq: ");
        for (int i = 0; i < logiSeq.size(); i++) {
            for (int j = 0; j < logiSeq.get(i).size(); j++) {
                System.out.print(logiSeq.get(i).get(j)+" -> ");
            }
            System.out.print(",   ");
        }
        System.out.println( );
    }

    private void printlogiSeq2(ArrayList<ArrayList<Integer>> logiSeq,int rv, int nv) {
        System.out.println("rv: ");
        for (int j = 0; j < logiSeq.get(rv).size(); j++) {
            System.out.print(logiSeq.get(rv).get(j)+" -> ");
        }
        System.out.println( );
        System.out.println("nv: ");
        for (int j = 0; j < logiSeq.get(nv).size(); j++) {
            System.out.print(logiSeq.get(nv).get(j)+" -> ");
        }
        System.out.println( );
    }

    private double getSurDual(double[] lambdaSur, double[] tadSur, double[] tssSur, double[] tafSur) {
        double obj=0;
        for (int i = 0; i < para.prodSet.size(); i++) {
            obj+=tafSur[i];
        }
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            obj+=lambdaSur[i]*(tadSur[i]-tssSur[i]);
        }
        return obj;
    }

    private double[] updateSurLambda(double[] gradSur, double[] lambdaSur, double stepsize) {
        double[] surLambda=new double[para.logiTaskSet.size()];
        for (int i = 0; i < surLambda.length; i++) {
            surLambda[i]=Math.max(0.0,lambdaSur[i]+stepsize*gradSur[i]);
        }
        return surLambda;
    }

    private double getSurNorm(double[] surGrad) {
        double norm=0;
        for (int i = 0; i < surGrad.length; i++) {
            norm+=surGrad[i]*surGrad[i];
        }
        return Math.sqrt(norm);
    }

    private double[] getSurGrad(double[] tadSur, double[] tssSur) {
        double[] surGrad=new double[para.logiTaskSet.size()];
        for (int i = 0; i < surGrad.length; i++) {
            surGrad[i]=tadSur[i]-tssSur[i];
        }
        return surGrad;
    }

    private double getSurDual() {
        double surDual = 0;
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            surDual += para.prodSet.get(prodSeq[p]).getTaf();
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    String key2 = prodSeq[p] + "_" + m + "_" + k;
                    double tss = para.prodSet.get(prodSeq[p]).getTaf()
                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
                            / para.lineSpeed;
                    double tad=PMKtoLogiTask.get(key2).getTad();
                    double lambda= PMKtoLogiTask.get(key2).getLambda();
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    surDual = surDual +( tad- tss) * lambda;
                }
            }
        }
        return surDual;
    }

//    public void run3_backup() throws ClassNotFoundException, IOException, IloException {
//        boolean flag = true;
//        bestLambdaForLowerBound = new double[logiTaskSet.size()];
//        double[] bestTssForLowerBound = new double[logiTaskSet.size()];
//        double[] bestTadForLowerBound = new double[logiTaskSet.size()];
//        // 为每一个logiTask设置其lambda（初始化拉格朗日乘子）
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
//                    PMKtoLogiTask.get(p + "_" + m + "_" + k).setLambda(Math.random()*0.);// (double)Math.random()*0.1
//                }
//            }
//        }
//        for (LogiTask logiTask : logiTaskSet) {
//            bestLambdaForLowerBound[logiTask.getId()] = 0.0;
//        }
//        iter = 0;
//        while (flag) {// (flag && iter<1000) {
//            iter++;
////			System.out.println("========iteration " + iter + ", incumbent bestRelaxObj="+bestRelaxObj+", bestObjUpBound="+bestObjUpBound +", bestIter="+bestIter+ " ============");
////			给定当前的lambda，求解lagrangian松弛问题 （最关键）
//            prodObj = solveProdProblem();// 子问题求解
//            logiObj = solveLogiProblem();// 子问题求解
//            LBNoLift = prodObj + logiObj;
//            double prodObjNoImprove=solveProdProblemNoLift();
//            LBLift =Math.max(LBLift,prodObjNoImprove+logiObj);
//            // 给定prodSeq得到两个算法下的可行目标值即上界
//            double objUpBound1 = LagrangHeuristic();
//            if (objUpBound1 < bestObjUpBound) {
//                bestObjUpBound = objUpBound1;
//            }
//            UBNoImprove=Math.min(objUpBound1,UBNoImprove);
//
//            double temp = LocalSearch(prodSeq,logiSeq);
//            if (temp < bestObjUpBound) {
//                bestObjUpBound = temp;
//                for (int ii = 0; ii < prodSeq.length; ii++) {
//                    bestProdSeq[ii] = prodSeq[ii];
//                }
//
//                bestLogiSeq.clear();
//                for (int ii = 0; ii < logiSeq.size(); ii++) {
//                    ArrayList<Integer> alist = new ArrayList<Integer>();
//                    for (int jj = 0; jj < logiSeq.get(ii).size(); jj++) {
//                        alist.add(logiSeq.get(ii).get(jj));
//                    }
//                    bestLogiSeq.add(alist);
//                }
//            }
//            if (LBNoLift > bestRelaxObj) {// 本次迭代松弛问题的目标值有改进(提升)
//                improved = true;
//                bestRelaxObj = LBNoLift;
//                for (LogiTask logiTask : logiTaskSet) {
//                    bestLambdaForLowerBound[logiTask.getId()] = logiTask.getLambda();
//                    bestTssForLowerBound[logiTask.getId()] = logiTask.getTss();
//                    bestTadForLowerBound[logiTask.getId()] = logiTask.getTad();
//                }
//                iter_no_improve = 0;
//            } else {
//                improved = false;
//                iter_no_improve++;
//            }
////          获得梯度
//            updateSubgradient();// update lambda with dualSolAndObj.solution
////          更新步长
//            updateStepsize();
//            // 输出信息
//            if (show) {
//                InfoOutput();
//            }
////			更新乘子
//            updateLambda();
//            if (iter_no_improve > max_iter_no_improve) {
//                flag = false;
//            }
//            if (iter > 1 && bestObjUpBound <= bestRelaxObj) {
//                bestRelaxObj = bestObjUpBound;
//                flag = false;
//            }
//            if (iter == 1 && bestRelaxObj - bestObjUpBound > 0.005) {
//                System.out.println("wrong");
//            }
////            System.out.println("iter="+iter+", LB="+bestRelaxObj+", UB="+bestObjUpBound);
//        }
////        System.out.println(tafLB[para.numProd-1]<para.prodSet.get(prodSeq[para.numProd-1]).getTaf());
////        System.out.println(para.prodSet.get(prodSeq[para.numProd-1]).getTaf());
////        System.out.println(" ");
//    }

    public void InfoOutput() throws IloException {
        System.out.println(">>>>>>>>>>>>>生产子问题信息<<<<<<<<<<<<");
        // 子问题信息输出
        System.out.println("production sequence:");
        for (int p = 0; p < para.numProd; p++) {
            System.out.print("-> " + prodSeq[p] + " ");
        }
        System.out.println();
        System.out.printf("%-16s\t", "ProdID");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10d", p);
        }
        System.out.println();

        System.out.printf("%-16s\t", "taf");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10.2f", para.prodSet.get(p).getTaf());
        }
        System.out.println();

//		System.out.printf("%-16s\t", "UB_Prod");
//		for (int p = 0; p < para.prodSet.size(); p++) {
//			System.out.printf("%10.2e", stepsizeUpperBound_Prod[p]);
//		}
//		System.out.println();

        System.out.println("prodObj=" + prodObj);

//		ProdProblem prodProblem = new ProdProblem(para);
//		prodProblem.buildModel();
//		prodProblem.getSolution();
//		System.out.println("optimal total completion time derived from MIP model: " + prodProblem.objVal);

        System.out.println(">>>>>>>>>>>>>物流子问题信息<<<<<<<<<<<<");
        // 子问题信息输出
        System.out.println("logistics sequence:");
//		for (AgvPath agvPath : optAgvPaths) {
//			System.out.print(agvPath);
//        }

        System.out.printf("%-16s\t", "LogiID");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10d", i);
        }
        System.out.println();

        System.out.printf("%-16s\t", "lambda");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.4f", para.ItoLogiTask.get(i).getLambda());
        }
        System.out.println();

        System.out.printf("%-16s\t", "tad");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.2f", para.ItoLogiTask.get(i).getTad());
        }
        System.out.println();

        System.out.printf("%-16s\t", "tss");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.2f", para.ItoLogiTask.get(i).getTss());
        }
        System.out.println();

        System.out.printf("%-16s\t", "subgrad");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.2f", para.ItoLogiTask.get(i).getSubgrad());
        }
        System.out.println();

//		System.out.printf("%-16s\t", "UB_Logi");
//		for (int i = 0; i < para.logiTaskSet.size(); i++) {
//			if (para.ItoLogiTask.get(i).getTad() < para.ItoLogiTask.get(i).getTss()) {// 该下标集合内各元素对应的lambda会减小
//				// 较小的步长以确保下一次的lambda不能为负数
//				System.out.printf("%10.4f", para.ItoLogiTask.get(i).getLambda()
//						/ (para.ItoLogiTask.get(i).getTss() - para.ItoLogiTask.get(i).getTad()));
//			} else {
//				System.out.printf("%10s", "Inf");
//			}
//		}
        System.out.println();
        System.out.println("logiObj=" + logiObj);
        // 松弛问题信息输出
        System.out.println(">>>>>>>>>>>>>松弛问题信息<<<<<<<<<<<<");
        System.out.println("stepsizeFactor=" + stepsizeFactor + ", iter=" + iter + ", bestRelaxObj=" + bestRelaxObj
                + ", bestObjUpBound=" + bestObjUpBound + ", stepsize=" + stepsize);
    }



//    private double solveProdProblem() throws IloException {
////		解松弛后的生产子问题，更新logiTask的tss
////		System.out.println("<<<<<<<<<<< SWPT rule For Production Problem >>>>>>>>>>>\n");
//
//        double[] sum_lambda_wrt_prod = new double[para.numProd];
//        double[] denominator_wrt_prod = new double[para.numProd];
//        double prodObj = 0;
//        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
//            int corresProdId = logiTask.getProd().getId();
//            sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + logiTask.getLambda();
//        }
//        for (int p = 0; p < para.numProd; p++) {
//            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
//                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
//        }
//        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
//        for (int p = 0; p < para.numProd; p++) {
//            ProdRatio prodRatio = new ProdRatio(p,
//                    (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
//            prodList.add(prodRatio);
//        }
//        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列
//
//
//        for (int p = 0; p < para.numProd; p++) {
//            prodSeq[p] = prodList.get(p).getId();
//        }
//
//        String key;
//        double[] taff=new double[para.prodSet.size()];
//        double timeCur =minTard;//Math.max(minTard,tardLB[prodSeq[0]]);// minTard;//para.moldChangeTime;//earliestReadyTime;// 时间游标
////        prodObj = 0;
//        double[] tempTaf = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
////			if (p>0) {
////				timeCur += tardLB[prodSeq[p]];
////			}
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
//                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
//                }
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    key = prodSeq[p] + "_" + m + "_" + k;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
//                        timeCur=Math.max( timeCur, tafLB[p] );
//                        para.prodSet.get(prodSeq[p]).setTaf(timeCur);
//                        taff[prodSeq[p]]=timeCur;
//                        tempTaf[prodSeq[p]] = timeCur;
//                    }
//                }
//            }
//            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
//                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
//        }
//
//
//        double[] latestTaf = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {
//            latestTaf[p] = tempTaf[p] + maxTaf - Arrays.stream(tempTaf).max().getAsDouble();//tempTaf[prodSeq[para.numProd - 1]];
//        }
//
//        for (int p = 0; p < para.numProd; p++) {
//            if (prodList.get(p).getRatio() < 0)//只有总权重小于0才会令其Taf尽可能大
//            {
//                para.prodSet.get(prodList.get(p).getId()).setTaf(latestTaf[prodList.get(p).getId()]);
//            }
//        }
//
//
//        // 基于新的taf更新tss
////        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
////            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
////                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
////                    double tss = para.prodSet.get(prodSeq[p]).getTaf()
////                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
////                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
////                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
////                            / para.lineSpeed;
//////                    PMKtoLogiTask.get(prodSeq[p] + "_" + m + "_" + k).setTss(tss);
////                }
////            }
////        }
//        // 算最终的目标函数
//        prodObj = 0;
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            prodObj = prodObj + para.prodSet.get(prodSeq[p]).getTaf();
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    String key2 = prodSeq[p] + "_" + m + "_" + k;
//                    double tss = para.prodSet.get(prodSeq[p]).getTaf()
//                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
//                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
//                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
//                            / para.lineSpeed;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    prodObj = prodObj - tss * PMKtoLogiTask.get(key2).getLambda();
//                }
//            }
//        }
//
//        return prodObj;
//    }

//    private double solveProdProblemNoLift(double[] lambdaSur) throws IloException {
////		解松弛后的生产子问题，更新logiTask的tss
////		System.out.println("<<<<<<<<<<< SWPT rule For Production Problem >>>>>>>>>>>\n");
//
//        double[] sum_lambda_wrt_prod = new double[para.numProd];
//        double[] denominator_wrt_prod = new double[para.numProd];
//        double prodObj = 0;
//        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
//            int corresProdId = logiTask.getProd().getId();
//            sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + lambdaSur[corresProdId];
//        }
//        for (int p = 0; p < para.numProd; p++) {
//            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
//                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
//        }
//        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
//        for (int p = 0; p < para.numProd; p++) {
//            ProdRatio prodRatio = new ProdRatio(p,
//                    (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
//            prodList.add(prodRatio);
//        }
//        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列
//
//
//        for (int p = 0; p < para.numProd; p++) {
//            prodSeq[p] = prodList.get(p).getId();
//        }
//
//        String key;
//        double timeCur = 0;//Math.max(minTard,tardLB[prodSeq[0]]);// minTard;//para.moldChangeTime;//earliestReadyTime;// 时间游标
////        prodObj = 0;
//        double[] tempTaf = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
////			if (p>0) {
////				timeCur += tardLB[prodSeq[p]];
////			}
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
//                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
//                }
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    key = prodSeq[p] + "_" + m + "_" + k;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
//                        tempTaf[prodSeq[p]] = timeCur;
//                    }
//                }
//            }
////            prodObj = prodObj + timeCur * para.objValFactor;
//            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
//                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
//        }
//
//
//        double[] latestTaf = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {
//            latestTaf[p] = tempTaf[p] + maxTaf - Arrays.stream(tempTaf).max().getAsDouble();//tempTaf[prodSeq[para.numProd - 1]];
//        }
//
//        for (int p = 0; p < para.numProd; p++) {
//            if (prodList.get(p).getRatio() < 0)//只有总权重小于0才会令其Taf尽可能大
//            {
//                taf_sur[prodList.get(p).getId()]=latestTaf[prodList.get(p).getId()];
//            }else {
//                taf_sur[prodList.get(p).getId()]=tempTaf[prodList.get(p).getId()];
//
//            }
//        }
//
//
//        // 基于新的taf更新tss
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
////            taf_sur[prodSeq[p]]=para.prodSet.get(prodSeq[p]).getTaf();
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    double tss = taf_sur[prodSeq[p]]
//                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
//                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
//                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
//                            / para.lineSpeed;
//                    tss_sur[PMKtoLogiTask.get(prodSeq[p] + "_" + m + "_" + k).getId()]=tss;
//                }
//            }
//        }
//        // 算最终的目标函数
//        prodObj = 0;
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            prodObj = prodObj + taf_sur[prodSeq[p]];
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    String key2 = prodSeq[p] + "_" + m + "_" + k;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    double tss=tss_sur[PMKtoLogiTask.get(key2).getId()];
//                    double lambda=lambdaSur[PMKtoLogiTask.get(key2).getId()];
//                    prodObj = prodObj - tss * lambda;
//                }
//            }
//        }
//
//        return prodObj;
//    }

//    private double solveProdProblemNoLift() throws IloException {
////		解松弛后的生产子问题，更新logiTask的tss
////		System.out.println("<<<<<<<<<<< SWPT rule For Production Problem >>>>>>>>>>>\n");
//
//        double[] sum_lambda_wrt_prod = new double[para.numProd];
//        double[] denominator_wrt_prod = new double[para.numProd];
//        double prodObj = 0;
//        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
//            int corresProdId = logiTask.getProd().getId();
//            sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + logiTask.getLambda();
//        }
//        for (int p = 0; p < para.numProd; p++) {
//            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
//                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
//        }
//        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
//        for (int p = 0; p < para.numProd; p++) {
//            ProdRatio prodRatio = new ProdRatio(p,
//                    (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
//            prodList.add(prodRatio);
//        }
//        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列
//
//
//        for (int p = 0; p < para.numProd; p++) {
//            prodSeq[p] = prodList.get(p).getId();
//        }
//
//        String key;
//        double timeCur = 0;//Math.max(minTard,tardLB[prodSeq[0]]);// minTard;//para.moldChangeTime;//earliestReadyTime;// 时间游标
////        prodObj = 0;
//        double[] tempTaf = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
////			if (p>0) {
////				timeCur += tardLB[prodSeq[p]];
////			}
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
//                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
//                }
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    key = prodSeq[p] + "_" + m + "_" + k;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
//                        para.prodSet.get(prodSeq[p]).setTaf(timeCur);
//                        tempTaf[prodSeq[p]] = timeCur;
//                    }
//                }
//            }
////            prodObj = prodObj + timeCur * para.objValFactor;
//            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
//                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
//        }
//
////        for (int p = 0; p < para.numProd; p++) {
////            double temp = para.prodSet.get(prodList.get(prodSeq[p]).getId()).getTaf();
////            if (temp < tafLB[p]) {
////                para.prodSet.get(prodList.get(prodSeq[p]).getId()).setTaf(tafLB[p]);
////            }
////            tempTaf[prodSeq[p]]=Math.max(temp,tafLB[p]);
////        }
//
//        double[] latestTaf = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {
//            latestTaf[p] = tempTaf[p] + maxTaf - Arrays.stream(tempTaf).max().getAsDouble();//tempTaf[prodSeq[para.numProd - 1]];
//        }
//
//        for (int p = 0; p < para.numProd; p++) {
//            if (prodList.get(p).getRatio() < 0)//只有总权重小于0才会令其Taf尽可能大
//            {
//                para.prodSet.get(prodList.get(p).getId()).setTaf(latestTaf[prodList.get(p).getId()]);
//            }
//        }
//
//
//        // 基于新的taf更新tss
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            taf_sur[prodSeq[p]]=para.prodSet.get(prodSeq[p]).getTaf();
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    double tss = para.prodSet.get(prodSeq[p]).getTaf()
//                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
//                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
//                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
//                            / para.lineSpeed;
//                    PMKtoLogiTask.get(prodSeq[p] + "_" + m + "_" + k).setTss(tss);
//                    tss_sur[PMKtoLogiTask.get(prodSeq[p] + "_" + m + "_" + k).getId()]=tss;
//                }
//            }
//        }
//        // 算最终的目标函数
//        prodObj = 0;
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            prodObj = prodObj + para.prodSet.get(prodSeq[p]).getTaf();
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    String key2 = prodSeq[p] + "_" + m + "_" + k;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    double tss=PMKtoLogiTask.get(key2).getTss();
//                    double lambda=PMKtoLogiTask.get(key2).getLambda();
//                    prodObj = prodObj - tss * lambda;
//                }
//            }
//        }
//
//        return prodObj;
//    }

    private double getLogiObj(double[] lambda, ArrayList<ArrayList<LogiTask>> logiSeq) {
        double obj=0;
        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                double tad=cur+ para.getDistance(para.depot, logiSeq.get(i).get(j).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                obj += tad * lambda[logiSeq.get(i).get(j).getId()];
                cur = tad;
                cur = cur + para.getDistance(para.depot, logiSeq.get(i).get(j).getMach())
                        / para.vehSpeed;
            }
        }
        return obj;
    }

    public ArrayList<ArrayList<LogiTask>> listLogiSeqCopy(ArrayList<ArrayList<LogiTask>> a){
        ArrayList<ArrayList<LogiTask>> b=new ArrayList<>();
        for (int i = 0; i < a.size(); i++) {
            b.add(new ArrayList<LogiTask>());
        }
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < a.get(i).size(); j++) {
                b.get(i).add(a.get(i).get(j));
            }
        }
        return b;
    }

//    private double solveLogiProblem_Ini(double[] lambda)  {
//
//        ArrayList<ArrayList<LogiTask>> cumb_logiSeq=new ArrayList<>();
//        for (int r = 0; r < para.numVeh; r++) {
//            cumb_logiSeq.add(new ArrayList<LogiTask>());
//        }
//
//        ArrayList<LogiTask> logiTaskList=new ArrayList<>();
//        for (LogiTask logiTask: para.logiTaskSet) {
//            logiTaskList.add(logiTask);
//        }
////        for (int i = 0; i < para.numVeh; i++) {
////            System.out.println(checkLambdaProcessTimeRatio1(logiTaskList));
////        }
//        Collections.sort(logiTaskList,new LogiTaskSortByLambdaProcessTimeRatio());
////        for (int i = 0; i < para.numVeh; i++) {
////            System.out.println(checkLambdaProcessTimeRatio1(logiTaskList));
////        }
//
//        double[] timeCur=new double[para.numVeh];
//        for (int i = 0; i < para.numVeh; i++) {
//            timeCur[i]=0;
//        }
//        int min_index=-1;
//        for (LogiTask logi : logiTaskList) {
//            double requiredDeliTime = logi.getChargingTime()
//                    + logi.getTransTime() / 2;
//            double min_val = Float.MAX_VALUE;
//            for (int r = 0; r < para.numVeh; r++) {
//                if (timeCur[r] + requiredDeliTime < min_val) {
//                    min_val = timeCur[r] + requiredDeliTime;
//                    min_index = r;
//                }
//            }
//            cumb_logiSeq.get(min_index).add(logi);
//            tad[logi.getId()] = timeCur[min_index] + logi.getChargingTime()
//                    + logi.getTransTime() / 2;
//            timeCur[min_index] = timeCur[min_index] + logi.getProcessTime();
//        }
//
//
//
//        ArrayList<ArrayList<LogiTask>> logiSeq = listLogiSeqCopy(cumb_logiSeq);
//
//
//        // 根据物流顺序返回tad的值
//        double[] tadd = new double[logiTaskSet.size()];
//        for (int i = 0; i < logiTaskSet.size(); i++) {
//            tadd[i] = Float.MAX_VALUE;
//        }
//        for (int i = 0; i < logiSeq.size(); i++) {
//            double baseTime = 0;
//            for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                LogiTask logiTask = logiSeq.get(i).get(j);
//                double temp = baseTime + logiTask.getChargingTime() + logiTask.getTransTime() / 2;
//                if (temp < tadd[logiTask.getId()]) {
//                    tadd[logiTask.getId()] = temp;
//                    tad_sur[logiTask.getId()] = temp;
//                }
//                baseTime = baseTime + logiTask.getProcessTime();
//            }
//        }
//        for (LogiTask logiTask : logiTaskSet) {
//            logiTask.setTad(tadd[logiTask.getId()]);
//        }
////
//        this.logiSeq.clear();
//        for (int r = 0; r < para.numVeh; r++) {
//            this.logiSeq.add(new ArrayList<Integer>());
//        }
//        for (int r = 0; r < logiSeq.size(); r++) {
//            for (int j=0;j<logiSeq.get(r).size();j++) {
//                this.logiSeq.get(r).add(logiSeq.get(r).get(j).getId());
//            }
//        }
//        return getLogiObj(lambda,logiSeq);
//
//    }

    private ArrayList<ArrayList<Integer>> getSurLogiSeq_Heur(double[] lambda, double targetObj,ArrayList<ArrayList<Integer>> logiSeq_ini)  {


        ArrayList<ArrayList<Integer>> bestLogiSeq=new ArrayList<>();


        ArrayList<LogiTask> logiTaskList=new ArrayList<>();
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            LogiTask logi=new LogiTask();
            logi.setId(i);
            logi.setLambda(lambda[i]);
            logi.setProcessingTime(para.ItoLogiTask.get(i).getProcessTime());
            logi.setTransTime(para.ItoLogiTask.get(i).getTransTime());
            logi.setChargingTime(para.ItoLogiTask.get(i).getChargingTime());
            logiTaskList.add(logi);
        }

        Collections.sort(logiTaskList, new LogiTaskSortByLambdaProcessTimeRatio());
//        for (int i = 0; i < para.numVeh; i++) {
//            System.out.println(checkLambdaProcessTimeRatio1(logiTaskList));
//        }

//        double[] timeCur=new double[para.numVeh];
//        for (int i = 0; i < para.numVeh; i++) {
//            timeCur[i]=0;
//        }
//        int min_index=-1;
//        for (int i=0;i<logiTaskList.size();i++) {
//            LogiTask logi=logiTaskList.get(i);
//            double requiredDeliTime = logi.getChargingTime()
//                    + logi.getTransTime() / 2;
//            double min_val = Float.MAX_VALUE;
//            for (int r = 0; r < para.numVeh; r++) {
//                if (timeCur[r] + requiredDeliTime < min_val) {
//                    min_val = timeCur[r] + requiredDeliTime;
//                    min_index = r;
//                }
//            }
//            logiSeq.get(min_index).add(logi.getId());
//            timeCur[min_index] = timeCur[min_index] + logi.getProcessTime();
//        }

        double bestObjVal=targetObj;
        boolean flag=true;
        Random rnd1 = new Random();
        Random rnd = new Random();
        Random rndOperator = new Random();
        int iter=0;
        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<>();

        while( iter<1){
//            System.out.println("first iteration "+iter);
            logiSeq.clear();
            for (int r = 0; r < para.numVeh; r++) {
                logiSeq.add(new ArrayList<Integer>());
            }
            int min_index=-1;
            for  (int i=0;i<logiTaskList.size();i++) {
                int logiTaskId=logiTaskList.get(i).getId();
                min_index = rnd1.nextInt(para.numVeh);
                logiSeq.get(min_index).add(logiTaskId);
            }

            double[] tad_temp=getSurTad(logiSeq);
            double temp=getSurLogiObj(lambda,tad_temp);
            if(temp<bestObjVal){
//                flag=true;
//                tad_sur=getTad(logiSeq);
//                this.logiSeq=getLogSeq(logiSeq);
                bestObjVal = temp;
            }

            int iter2 = 0;

            while (iter2 < 1) {
//                System.out.println("second iteration "+iter2);
                int ov;
                int op;
                while (true) {
                    int a = rnd.nextInt(para.numVeh);
                    if (logiSeq.get(a).size() >= 1) {
                        ov = a;
                        op = rnd.nextInt(logiSeq.get(a).size());
                        break;
                    }
                }


                boolean improved = false;
                int iter_inner=0;
                while (!improved && iter_inner < 5) {
//                    System.out.println("third iteration "+iter_inner);
                    switch (rndOperator.nextInt(2)) {
                        case 0: {//relocate
//                            System.out.println("case0: relocate");
                            int nv;
//                            printlogiSeq(logiSeq);
//                            System.out.println("ov="+ov+", op="+op);
                            while (true){
                                nv=rnd.nextInt(para.numVeh);
                                if(nv!=ov){
//                                    System.out.println("relocate");
                                    break;
                                }
                            }
//                            System.out.println("nv="+nv);
//                            System.out.println("before change:");
//                            printlogiSeq(logiSeq);
//                            printRatio(formFullLogiSeq(logiSeq,lambda));
//                            checkLambdaProcessTimeRatio(formFullLogiSeq(logiSeq,lambda));
                            ArrayList<ArrayList<Integer>> neigh_logiSeq=logiRelocateOperator(lambda,logiSeq, ov, op, nv);
//                            System.out.println("after change:");
//                            printlogiSeq(neigh_logiSeq);
//                            printRatio(formFullLogiSeq(neigh_logiSeq,lambda));
//                            checkLambdaProcessTimeRatio(formFullLogiSeq(neigh_logiSeq,lambda));

                            tad_temp=getSurTad(neigh_logiSeq);
                            double testObjVal = getSurLogiObj(lambda,tad_temp);
                            if (testObjVal < bestObjVal) {
                                bestObjVal = testObjVal;
                                logiSeq = listCopy2(neigh_logiSeq);
//                                printlogiSeq(logiSeq);
                                improved = true;
                            }
                            break;
                        }
                        case 1: {//swap
                            if(!checkSwapFeasible(ov,logiSeq)){
                                break;
                            }
//                            System.out.println("case1: swap");
                            int nv,np;
//                            printlogiSeq(logiSeq);
//                            System.out.println("ov="+ov+", op="+op);
                            while (true){
                                nv=rnd.nextInt(para.numVeh);
                                if(nv!=ov && logiSeq.get((nv)).size()>=1){
                                    np = rnd.nextInt(logiSeq.get((nv)).size());
//                                    System.out.println("swap");
                                    break;
                                }
                            }
//                            System.out.println("nv="+nv+", np="+np);
//                            System.out.println("before change:");
//                            printlogiSeq(logiSeq);
//                            printRatio(formFullLogiSeq(logiSeq,lambda));
                            ArrayList<ArrayList<Integer>> neigh_logiSeq=logiSwapOperator(lambda,logiSeq, ov, op, nv, np);
//                            System.out.println("after change:");
//                            printlogiSeq(neigh_logiSeq);
//                            printRatio(formFullLogiSeq(neigh_logiSeq,lambda));
//                            checkLambdaProcessTimeRatio(formFullLogiSeq(neigh_logiSeq,lambda));

                            tad_temp=getSurTad(neigh_logiSeq);
                            double testObjVal = getSurLogiObj(lambda,tad_temp);
                            if (testObjVal < bestObjVal) {
                                bestObjVal = testObjVal;
                                logiSeq = listCopy2(neigh_logiSeq);
//                                printlogiSeq(logiSeq);
                                improved = true;
                            }
                            break;
                        }
                    }
                    iter_inner++;
                }
                iter2++;
            }
            iter++;
        }
        if (bestObjVal<targetObj){
            bestLogiSeq=listCopy2(logiSeq);
//            System.out.println("logistics subproblem solved: find a solution better than targetObj");
        }else{
            bestLogiSeq=listCopy2(logiSeq_ini);
//            System.out.println("logistics subproblem solved: not find better solution");
        }
//        System.out.println("found better logiSeq: "+ flag);
        return bestLogiSeq;

    }

    private boolean checkSwapFeasible(int ov, ArrayList<ArrayList<Integer>> logiSeq) {
        boolean flag=false;
        for (int i = 0; i < para.numVeh; i++) {
            if(i!=ov){
                if(logiSeq.get(i).size()>0){
                    flag=true;
                }
            }
        }
        return flag;
    }
//    private double solveLogiProblem_Heur(double[] lambda, double targetObj)  {
//        boolean flag=false;
//        double bestObjVal=targetObj;
//
//        ArrayList<ArrayList<LogiTask>> cumb_logiSeq=new ArrayList<>();
//        for (int r = 0; r < para.numVeh; r++) {
//            cumb_logiSeq.add(new ArrayList<LogiTask>());
//        }
//
//        ArrayList<LogiTask> logiTaskList=new ArrayList<>();
//        for (LogiTask logiTask: para.logiTaskSet) {
//            logiTaskList.add(logiTask);
//        }
//        Collections.sort(logiTaskList,new LogiTaskSortByLambdaProcessTimeRatio());
//        ArrayList<ArrayList<LogiTask>> logiSeq=new ArrayList<>();
//
//        for (int ii = 0; ii < 1; ii++) {
//
//            Random rnd1 = new Random();
//            int min_index;
//            for (LogiTask logi : logiTaskList) {
//                min_index = rnd1.nextInt(para.numVeh);
//                cumb_logiSeq.get(min_index).add(logi);
//            }
//
//            logiSeq = listLogiSeqCopy(cumb_logiSeq);
//            double temp=getLogiObj(lambda,logiSeq);
//            if(temp<bestObjVal){
//                flag=true;
//                tad_sur=getTad(logiSeq);
//                this.logiSeq=getLogSeq(logiSeq);
//                bestObjVal = temp;
//            }
//
//            int iter = 0;
//            Random rnd = new Random();
//            while (iter < 100) {
//                int ov;
//                int op;
//                while (true) {
//                    int a = rnd.nextInt(para.numVeh);
//                    if (logiSeq.get(a).size() >= 1) {
//                        ov = a;
//                        op = rnd.nextInt(logiSeq.get(a).size());
//                        break;
//                    }
//                }
//
//                Random rndOperator = new Random();
//                boolean improved = false;
//                int iter_inner=0;
//                while (!improved && iter_inner < 2000) {
//                    iter_inner++;
//                    switch (rndOperator.nextInt(2)) {
//                        case 0: {//relocate
//                            int nv;
//                            while (true){
//                                nv=rnd.nextInt(para.numVeh);
//                                if(nv!=ov){
//                                    break;
//                                }
//                            }
//                            ArrayList<ArrayList<LogiTask>> neigh_logiSeq = listLogiSeqCopy(logiSeq);
//                            logiRelocateOperator(neigh_logiSeq, ov, op, nv);
//
//
//                            double testObjVal = getLogiObj(lambda,neigh_logiSeq);
//                            if (testObjVal < bestObjVal) {
//                                flag=true;
//                                tad_sur=getTad(neigh_logiSeq);
//                                this.logiSeq=getLogSeq(neigh_logiSeq);
//                                bestObjVal = testObjVal;
//                                logiSeq = neigh_logiSeq;
//                                improved = true;
//                            }
//                            break;
//                        }
//                        case 1: {//swap
//                            int nv,np;
//                            while (true){
//                                nv=rnd.nextInt(para.numVeh);
//                                if(nv!=ov && logiSeq.get((nv)).size()>=1){
//                                    np = rnd.nextInt(logiSeq.get((nv)).size());
//                                    break;
//                                }
//                            }
//                            ArrayList<ArrayList<LogiTask>> neigh_logiSeq = listLogiSeqCopy(logiSeq);
//                            logiSwapOperator(neigh_logiSeq, ov, op, nv, np);
//
//                            double testObjVal = getLogiObj(lambda,neigh_logiSeq);
//                            if (testObjVal < bestObjVal) {
//                                flag=true;
//                                tad_sur=getTad(neigh_logiSeq);
//                                this.logiSeq=getLogSeq(neigh_logiSeq);
//                                bestObjVal = testObjVal;
//                                logiSeq = neigh_logiSeq;
//                                improved = true;
//                            }
//                            break;
//                        }
//                    }
//                }
//                iter++;
//            }
//        }
//        System.out.println("found better logiSeq: "+ flag);
//        return bestObjVal;
//
//    }

    private ArrayList<ArrayList<Integer>> getLogSeq(ArrayList<ArrayList<LogiTask>> neighLogiSeq) {
        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<>();
        for (int i = 0; i < neighLogiSeq.size(); i++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < neighLogiSeq.size(); i++) {
            for (int j = 0; j < neighLogiSeq.get(i).size(); j++) {
                logiSeq.get(i).add(neighLogiSeq.get(i).get(j).getId());
            }
        }

        return logiSeq;
    }

    private double[] getTad(ArrayList<ArrayList<LogiTask>> logiSeq) {
        double tad[]=new double[para.logiTaskSet.size()];
        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                tad[logiSeq.get(i).get(j).getId()]=cur+ para.getDistance(para.depot, logiSeq.get(i).get(j).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = tad[logiSeq.get(i).get(j).getId()] + para.getDistance(para.depot, logiSeq.get(i).get(j).getMach())
                        / para.vehSpeed;
            }
        }
        return tad;
    }

//    private double solveLogiProblem_Heur(double targetObj)  {
//
//        double bestObjVal=Float.MAX_VALUE;
//
//        ArrayList<ArrayList<LogiTask>> cumb_logiSeq=new ArrayList<>();
//        for (int r = 0; r < para.numVeh; r++) {
//            cumb_logiSeq.add(new ArrayList<LogiTask>());
//        }
//
//        ArrayList<LogiTask> logiTaskList=new ArrayList<>();
//        for (LogiTask logiTask: para.logiTaskSet) {
//            logiTaskList.add(logiTask);
//        }
////        for (int i = 0; i < para.numVeh; i++) {
////            System.out.println(checkLambdaProcessTimeRatio1(logiTaskList));
////        }
//        Collections.sort(logiTaskList,new LogiTaskSortByLambdaProcessTimeRatio());
////        for (int i = 0; i < para.numVeh; i++) {
////            System.out.println(checkLambdaProcessTimeRatio1(logiTaskList));
////        }
//        ArrayList<ArrayList<LogiTask>> logiSeq=new ArrayList<>();
//
//        for (int ii = 0; ii < 1; ii++) {
//
//            Random rnd1 = new Random();
//            int min_index;
//            for (LogiTask logi : logiTaskList) {
//                min_index = rnd1.nextInt(para.numVeh);
//                cumb_logiSeq.get(min_index).add(logi);
//            }
//
//            logiSeq = listLogiSeqCopy(cumb_logiSeq);
//            bestObjVal = Math.min(getLogiObj(logiSeq), bestObjVal);
//
//            int iter = 0;
//            Random rnd = new Random();
//            while (iter < 2) {
//                int ov;
//                int op;
//                while (true) {
//                    int a = rnd.nextInt(para.numVeh);
//                    if (logiSeq.get(a).size() >= 1) {
//                        ov = a;
//                        op = rnd.nextInt(logiSeq.get(a).size());
//                        break;
//                    }
//                }
//
//                Random rndOperator = new Random();
//                boolean improved = false;
//                int iter_inner=0;
//                while (!improved && iter_inner < 20) {
//                    iter_inner++;
//                    switch (rndOperator.nextInt(2)) {
//                        case 0: {//relocate
//                            int nv;
//                            while (true){
//                                nv=rnd.nextInt(para.numVeh);
//                                if(nv!=ov){
//                                    break;
//                                }
//                            }
//                            ArrayList<ArrayList<LogiTask>> neigh_logiSeq = listLogiSeqCopy(logiSeq);
////                            if (!checkLambdaProcessTimeRatio(neigh_logiSeq)) {
////                                System.out.println("wrong");
////                            }
//                            logiRelocateOperator(neigh_logiSeq, ov, op, nv);
////                            if (!checkLambdaProcessTimeRatio(neigh_logiSeq)) {
////                                System.out.println("wrong");
////                            }
////                            checkLambdaProcessTimeRatio(neigh_logiSeq);
//                            double testObjVal = getLogiObj(neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq = neigh_logiSeq;
//                                improved = true;
//                            }
//                            break;
//                        }
//                        case 1: {//swap
//                            int nv,np;
//                            while (true){
//                                nv=rnd.nextInt(para.numVeh);
//                                if(nv!=ov && logiSeq.get((nv)).size()>=1){
//                                    np = rnd.nextInt(logiSeq.get((nv)).size());
//                                    break;
//                                }
//                            }
//                            ArrayList<ArrayList<LogiTask>> neigh_logiSeq = listLogiSeqCopy(logiSeq);
////                            if (!checkLambdaProcessTimeRatio(neigh_logiSeq)) {
////                                System.out.println("wrong");
////                            }
//                            logiSwapOperator(neigh_logiSeq, ov, op, nv, np);
////                            if (!checkLambdaProcessTimeRatio(neigh_logiSeq)) {
////                                System.out.println("wrong");
////                            }
////                            checkLambdaProcessTimeRatio(neigh_logiSeq);
//                            double testObjVal = getLogiObj(neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq = neigh_logiSeq;
//                                improved = true;
//                            }
//                            break;
//                        }
//                    }
//                }
//                iter++;
//            }
//        }
//        if(bestObjVal<targetObj) {
//
//            // 根据物流顺序返回tad的值
//            double[] tadd = new double[logiTaskSet.size()];
//            for (int i = 0; i < logiTaskSet.size(); i++) {
//                tadd[i] = Float.MAX_VALUE;
//            }
//            for (int i = 0; i < logiSeq.size(); i++) {
//                double baseTime = 0;
//                for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                    LogiTask logiTask = logiSeq.get(i).get(j);
//                    double temp = baseTime + logiTask.getChargingTime() + logiTask.getTransTime() / 2;
//                    if (temp < tadd[logiTask.getId()]) {
//                        tadd[logiTask.getId()] = temp;
//                    }
//                    baseTime = baseTime + logiTask.getProcessTime();
//                }
//            }
//            for (LogiTask logiTask : logiTaskSet) {
//                logiTask.setTad(tadd[logiTask.getId()]);
//            }
////
//            this.logiSeq.clear();
//            for (int r = 0; r < para.numVeh; r++) {
//                this.logiSeq.add(new ArrayList<Integer>());
//            }
//            for (int r = 0; r < logiSeq.size(); r++) {
//                for (int j=0;j<logiSeq.get(r).size();j++) {
//                    this.logiSeq.get(r).add(logiSeq.get(r).get(j).getId());
//                }
//            }
//        }
//        return Math.min(bestObjVal,targetObj);
//
//    }

    public ArrayList<ArrayList<Integer>> logiSwapOperator(double[] lambda, ArrayList<ArrayList<Integer>> input_logiSeq, int ov, int op, int nv, int np){
        ArrayList<ArrayList<LogiTask>> new_logiSeq=formFullLogiSeq(input_logiSeq, lambda);

//        for (int i = 0; i < para.numVeh; i++) {
//            checkLambdaProcessTimeRatio1(new_logiSeq.get(i));
//        }

        LogiTask logi1=new_logiSeq.get(ov).get(op);
        LogiTask logi2=new_logiSeq.get(nv).get(np);
        new_logiSeq.get(ov).remove(op);
        new_logiSeq.get(nv).remove(np);

//        for (int i = 0; i < para.numVeh; i++) {
//            checkLambdaProcessTimeRatio1(new_logiSeq.get(i));
//        }

        new_logiSeq.get(nv).add(logi1);
        new_logiSeq.get(ov).add(logi2);


        Collections.sort(new_logiSeq.get(nv),new LogiTaskSortByLambdaProcessTimeRatio());
        Collections.sort(new_logiSeq.get(ov),new LogiTaskSortByLambdaProcessTimeRatio());

//        for (int i = 0; i < para.numVeh; i++) {
//            checkLambdaProcessTimeRatio1(new_logiSeq.get(i));
//        }

        ArrayList<ArrayList<Integer>> new_logiSeq2 = logiSeqImprove(new_logiSeq,lambda);

        ArrayList<ArrayList<Integer>> logiSeq=listCopy2(new_logiSeq2);
        return logiSeq;
    }


    public ArrayList<ArrayList<Integer>> logiRelocateOperator(double[] lambda, ArrayList<ArrayList<Integer>> input_logiSeq, int ov, int op, int nv){
        ArrayList<ArrayList<LogiTask>> new_logiSeq=formFullLogiSeq(input_logiSeq, lambda);

//        for (int i = 0; i < para.numVeh; i++) {
//            checkLambdaProcessTimeRatio1(new_logiSeq.get(i));
//        }

        LogiTask logi=new_logiSeq.get(ov).get(op);
        new_logiSeq.get(ov).remove(op);
//        for (int i = 0; i < para.numVeh; i++) {
//            checkLambdaProcessTimeRatio1(new_logiSeq.get(i));
//        }
        new_logiSeq.get(nv).add(logi);

        Collections.sort(new_logiSeq.get(nv),new LogiTaskSortByLambdaProcessTimeRatio());
//        for (int i = 0; i < para.numVeh; i++) {
//            checkLambdaProcessTimeRatio1(new_logiSeq.get(i));
//        }

        ArrayList<ArrayList<Integer>> new_logiSeq2 = logiSeqImprove(new_logiSeq,lambda);

        ArrayList<ArrayList<Integer>> logiSeq=listCopy2(new_logiSeq2);

        return logiSeq;
    }

    private ArrayList<ArrayList<Integer>> logiSeqImprove(ArrayList<ArrayList<LogiTask>> newLogiSeq, double[] lambda) {

        ArrayList<ArrayList<Integer>> logiSeq=extractLogiTaskID(newLogiSeq);
//        int flag;

//        System.out.println("before impprove:");
//        printlogiSeq(logiSeq);
//        printArray("ts ",ts);
//        printArray("tc ",tc);
//        printRatio(formFullLogiSeq(logiSeq,lambda));

//        flag=findStartTimeAndCompletionTime(tc,ts,lastId);
        boolean flag=true;
//        Random rnd=new Random();
        while(flag){
//            int rv=rnd.nextInt(para.numVeh);
//            while (rv==flag){
//                rv=rnd.nextInt(para.numVeh);
            double[] tc = getSurCompletionTimeForParallelMachineSchedule(logiSeq);
            double[] ts = getSurStartTimeForParallelMachineSchedule(logiSeq);
            int[] lastId=getSurLastLogiTaskId(logiSeq);
//            }
            int ov=-1;
            int nv=-1;
            double minTc=Float.MAX_VALUE;
            double maxTs=Float.MIN_VALUE;
            int minTcVeh=-1;
            int maxTsVeh=-1;

            for (int i = 0; i < lastId.length; i++) {
                if(lastId[i]!=-1) {
                    if (tc[lastId[i]] < minTc) {
                        minTc = tc[lastId[i]];
                        minTcVeh=i;
                    }
                    if (ts[lastId[i]] > maxTs) {
                        maxTs = ts[lastId[i]];
                        maxTsVeh=i;
                    }
                }
            }

//            System.out.println("min"+minTcVeh);
//            System.out.println("max"+maxTsVeh);
//            System.out.println("lastmin"+lastId[minTcVeh]);
//            System.out.println("lastmax"+lastId[maxTsVeh]);
            if(tc[lastId[minTcVeh]]<ts[lastId[maxTsVeh]]){
                ov=maxTsVeh;
                nv=minTcVeh;
                for (int i = 0; i < lastId.length; i++) {
                    if(lastId[i]==-1) {
                        nv=i;
                    }
                }

    //            printlogiSeq2(logiSeq,rv,flag);
    //            System.out.println("rv lastID "+logiSeq.get(rv).getLast());
    //            System.out.println("rv lastID start time"+ts[logiSeq.get(rv).getLast()]);
    //            System.out.println("rv lastID comp time"+tc[logiSeq.get(rv).getLast()]);
    //            System.out.println(flag);
    //            System.out.println("nv lastID "+logiSeq.get(flag).getLast());
    //            System.out.println("nv lastID start time"+ts[logiSeq.get(flag).getLast()]);
    //            System.out.println("nv lastID comp time"+tc[logiSeq.get(flag).getLast()]);

                int a=logiSeq.get(ov).getLast();
                logiSeq.get(ov).removeLast();
                logiSeq.get(nv).add(a);
                ArrayList<ArrayList<LogiTask>> new_logiSeq=formFullLogiSeq(logiSeq,lambda);
                Collections.sort(new_logiSeq.get(nv),new LogiTaskSortByLambdaProcessTimeRatio());
                logiSeq=extractLogiTaskID(new_logiSeq);
    //            tc = getSurCompletionTimeForParallelMachineSchedule(logiSeq);
    //            ts = getSurStartTimeForParallelMachineSchedule(logiSeq);
    //            lastId=getSurLastLogiTaskId(logiSeq);
    //            flag=findStartTimeAndCompletionTime(tc,ts,lastId);
            }else {
                flag=false;
            }


        }

//        System.out.println("after impprove:");
//        printlogiSeq(logiSeq);
//        printArray("ts ",ts);
//        printArray("tc ",tc);
//        printRatio(formFullLogiSeq(logiSeq,lambda));

        return logiSeq;
    }

    private int findStartTimeAndCompletionTime(double[] tc, double[] ts, int[] lastId) {
        ArrayList<Integer> a=new ArrayList<>();
        for (int i = 0; i < para.numVeh; i++) {
            for (int j = 0; j < para.logiTaskSet.size(); j++) {
                if(lastId[i]==-1){
                    return i;
                } else if(tc[lastId[i]]<ts[j]){
                    a.add(i);
                }
            }
        }
        int g=-1;
        if(a.size()>0){
            double k=Float.MAX_VALUE;
            for (int i = 0; i < a.size(); i++) {
                if(tc[lastId[a.get(i)]]<k){
                    k=tc[lastId[a.get(i)]];
                    g=a.get(i);
                }
            }
        }else {
            g=-1;
        }

        return g;
    }

    private int[] getSurLastLogiTaskId(ArrayList<ArrayList<Integer>> logiSeq) {
        int[] last= new int[para.numVeh];
        for (int i = 0; i < logiSeq.size(); i++) {
            if(logiSeq.get(i).size()==0){
                last[i]=-1;
            }else {
                last[i]=logiSeq.get(i).getLast();
            }
        }

        return last;
    }

    private double[] getSurCompletionTimeForParallelMachineSchedule(ArrayList<ArrayList<Integer>> logiSeq) {

        double[] tc=new double[para.logiTaskSet.size()];

        for (int i = 0; i < logiSeq.size(); i++) {
            double baseTime = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {
                int logiTaskId = logiSeq.get(i).get(j);
                tc[logiTaskId]=baseTime + para.ItoLogiTask.get(logiTaskId).getProcessTime() ;
                baseTime = tc[logiTaskId];
            }
        }

        return tc;
    }

    private double[] getSurStartTimeForParallelMachineSchedule(ArrayList<ArrayList<Integer>> logiSeq) {

        double[] ts=new double[para.logiTaskSet.size()];

        for (int i = 0; i < logiSeq.size(); i++) {
            double baseTime = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {
                int logiTaskId = logiSeq.get(i).get(j);
                ts[logiTaskId] = baseTime;
                baseTime += para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }
        }

        return ts;
    }

    private void printRatio(ArrayList<ArrayList<LogiTask>> newLogiSeq) {
        for (int i = 0; i < newLogiSeq.size(); i++) {
            System.out.print("veh "+i+":");
            for (int j = 0; j < newLogiSeq.get(i).size(); j++) {
                System.out.printf("%14.6f",newLogiSeq.get(i).get(j).getLambda()/newLogiSeq.get(i).get(j).getProcessTime());
                System.out.print("   ");
            }
            System.out.println();
        }
    }

    private ArrayList<ArrayList<Integer>> extractLogiTaskID(ArrayList<ArrayList<LogiTask>> inputLogiSeq) {
        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<>();
        for (int i = 0; i < para.numVeh; i++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < inputLogiSeq.size(); i++) {
            for (int j = 0; j < inputLogiSeq.get(i).size(); j++) {
                logiSeq.get(i).add(inputLogiSeq.get(i).get(j).getId());
            }
        }
        return logiSeq;
    }

    private ArrayList<ArrayList<LogiTask>> formFullLogiSeq(ArrayList<ArrayList<Integer>> inputLogiSeq, double[] lambda) {
        ArrayList<ArrayList<LogiTask>> logiSeq=new ArrayList<>();
        for (int i = 0; i < para.numVeh; i++) {
            logiSeq.add(new ArrayList<LogiTask>());
        }
        for (int i = 0; i < inputLogiSeq.size(); i++) {
            for (int j = 0; j < inputLogiSeq.get(i).size(); j++) {
                int id=inputLogiSeq.get(i).get(j);

                LogiTask logi=new LogiTask();
                logi.setId(inputLogiSeq.get(i).get(j));
                logi.setLambda(lambda[id]);
                logi.setProcessingTime(para.ItoLogiTask.get(id).getProcessTime());
                logi.setTransTime(para.ItoLogiTask.get(id).getTransTime());
                logi.setChargingTime(para.ItoLogiTask.get(id).getChargingTime());
                logiSeq.get(i).add(logi);
            }
        }

        return logiSeq;
    }

//    private double solveLogiProblem_Heur2(double targetObj)  {
////        given lambda,solve LR-MH
////		启发式解物流子问题，更新logiTask的tad
////        ColGen colGen = new ColGen(para);
////        colGen.run();
////        double logiObj = colGen.getObj();// 列生成获得的obj是不带 sum_i (lambda_i*TransTime_(0,i))的
////
////        optAgvPaths = colGen.rstMasterProblem.getOptAgvPaths();// ArrayList<AgvPath>类型
//
////        prodSeq
//        double bestObjVal=Float.MAX_VALUE;
//        ArrayList<ArrayList<Integer>> cumb_logiSeq=new ArrayList<ArrayList<Integer>>();
//        for (int r = 0; r < para.numVeh; r++) {
//            cumb_logiSeq.add(new ArrayList<Integer>());
//        }
//
//        ArrayList<LogiTaskIdLambda> logiTaskIdLambdaList=new ArrayList<>();
//        for (LogiTask logiTask: para.logiTaskSet) {
//            LogiTaskIdLambda logiTaskIdLambda=new LogiTaskIdLambda(logiTask.getId(),logiTask.getLambda());
//            logiTaskIdLambdaList.add(logiTaskIdLambda);
//        }
//        Collections.sort(logiTaskIdLambdaList,new LogiTaskSortByLambda());
//        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<>();
//        for (int ii = 0; ii < 100; ii++) {
//
//            double[] timeCur = new double[para.numVeh];
//            double[] tad = new double[para.logiTaskSet.size()];
//            for (int r = 0; r < cumb_logiSeq.size(); r++) {
//                timeCur[r] = 0;
//            }
//            int min_index = 0;
//            for (LogiTaskIdLambda logi : logiTaskIdLambdaList) {
//                int logiTaskId = logi.id;
//                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
//                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
//                Random rnd = new Random();
////            double min_val = Float.MAX_VALUE;
////            for (int r = 0; r < para.numVeh; r++) {
////                if (timeCur[r] + requiredDeliTime < min_val) {
////                    min_val = timeCur[r] + requiredDeliTime;
////                    min_index = r;
////                }
////            }
//                min_index = rnd.nextInt(para.numVeh);
//                cumb_logiSeq.get(min_index).add(logiTaskId);
//                tad[logiTaskId] = timeCur[min_index] + requiredDeliTime;
//                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
//            }
//
//            logiSeq = listCopy2(cumb_logiSeq);
//            bestObjVal = Math.min(getLogiObj(logiSeq), bestObjVal);
//
//            int iter = 0;
//            Random rnd = new Random();
//            while (iter < 2000 && bestObjVal >= targetObj) {
////                System.out.println('j');
//                int ov;
//                int op;
//                while (true) {
//                    int a = rnd.nextInt(para.numVeh);
//                    if (logiSeq.get(a).size() >= 1) {
//                        ov = a;
//                        op = rnd.nextInt(logiSeq.get(a).size());
//                        break;
//                    }
//                }
//
//
//                int o_LogiTaskId = logiSeq.get(ov).get(op);
//                double o_Lambda = para.ItoLogiTask.get(o_LogiTaskId).getLambda();
//                double obefore_Lambda, oafter_Lambda;
//                if (logiSeq.get(ov).size() == 1) {
//                    obefore_Lambda = Float.MIN_VALUE;
//                    oafter_Lambda = Float.MAX_VALUE;
//                } else {
//                    if (op == 0) {
//                        obefore_Lambda = Float.MIN_VALUE;
//                        int oafter_LogiTaskId = logiSeq.get(ov).get(op + 1);
//                        oafter_Lambda = para.ItoLogiTask.get(oafter_LogiTaskId).getLambda();
//                    } else if (op == logiSeq.get(ov).size() - 1) {
//                        oafter_Lambda = Float.MAX_VALUE;
//                        int obefore_LogiTaskId = logiSeq.get(ov).get(op - 1);
//                        obefore_Lambda = para.ItoLogiTask.get(obefore_LogiTaskId).getLambda();
//                    } else {
//                        int obefore_LogiTaskId = logiSeq.get(ov).get(op - 1);
//                        int oafter_LogiTaskId = logiSeq.get(ov).get(op + 1);
//                        obefore_Lambda = para.ItoLogiTask.get(obefore_LogiTaskId).getLambda();
//                        oafter_Lambda = para.ItoLogiTask.get(oafter_LogiTaskId).getLambda();
//                    }
//                }
//
//                HashSet<int[]> feasBeforeInsertPos = new HashSet<>();//插入在指定位置之后
//                HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
//                HashSet<int[]> feasNewInsertPos = new HashSet<>();//插入在队列末尾
//                HashSet<int[]> feasSwapPos = new HashSet<>();//插入在队列末尾
//
//                for (int i = 0; i < logiSeq.size(); i++) {
//                    if (i != ov) {
//                        if (logiSeq.get(i).size() == 0) {
//                            int[] fea = {i, 0};
//                            feasNewInsertPos.add(fea);
//                        } else {
//                            for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                                int n_LogiTaskId = logiSeq.get(i).get(j);
//                                double n_Lambda = para.ItoLogiTask.get(n_LogiTaskId).getLambda();
//                                double nbefore_Lambda = 0, nafter_Lambda = 0;
//                                if (n_LogiTaskId != o_LogiTaskId) {
//                                    if (logiSeq.get(i).size() == 1) {
//                                        nbefore_Lambda = Float.MIN_VALUE;
//                                        nafter_Lambda = Float.MAX_VALUE;
//                                    } else {
//                                        if (j == 0) {
//                                            nbefore_Lambda = Float.MIN_VALUE;
//                                            int nafter_LogiTaskId = logiSeq.get(i).get(j + 1);
//                                            nafter_Lambda = para.ItoLogiTask.get(nafter_LogiTaskId).getLambda();
//                                        } else if (j == logiSeq.get(i).size() - 1) {
//                                            nafter_Lambda = Float.MAX_VALUE;
//                                            int nbefore_LogiTaskId = logiSeq.get(i).get(j - 1);
//                                            nbefore_Lambda = para.ItoLogiTask.get(nbefore_LogiTaskId).getLambda();
//                                        } else {
//                                            int nbefore_LogiTaskId = logiSeq.get(i).get(j - 1);
//                                            int nafter_LogiTaskId = logiSeq.get(i).get(j + 1);
//                                            nbefore_Lambda = para.ItoLogiTask.get(nbefore_LogiTaskId).getLambda();
//                                            nafter_Lambda = para.ItoLogiTask.get(nafter_LogiTaskId).getLambda();
//                                        }
//                                    }
//                                    int[] fea = {i, j};
//                                    if (o_Lambda <= n_Lambda && o_Lambda >= nbefore_Lambda) {
//                                        feasBeforeInsertPos.add(fea);
//                                    }
//                                    if (j == logiSeq.get(i).size() - 1) {
//                                        if (o_Lambda > n_Lambda) {
//                                            feasEndInsertPos.add(fea);
//                                        }
//                                    }
//                                    if (o_Lambda >= nbefore_Lambda && o_Lambda <= nafter_Lambda && n_Lambda >= obefore_Lambda && n_Lambda <= oafter_Lambda) {
//                                        feasSwapPos.add(fea);
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//
//
//                Random rndOperator = new Random();
//                boolean improved = false;
//                while (!improved && (!feasSwapPos.isEmpty() || !feasBeforeInsertPos.isEmpty() || !feasEndInsertPos.isEmpty() || !feasNewInsertPos.isEmpty())) {
//                    switch (rndOperator.nextInt(4)) {
//                        case 0:
//                            if (!feasBeforeInsertPos.isEmpty()) {
//                                Iterator<int[]> iterator = feasBeforeInsertPos.iterator();
//                                int[] fea = iterator.next();
//                                feasBeforeInsertPos.remove(fea);
//                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                                beforeInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                                if (!checkLambda(neigh_logiSeq)) {
//                                    System.out.println("wrong");
//                                }
//                                double testObjVal = getLogiObj(neigh_logiSeq);
//                                if (testObjVal <= bestObjVal) {
//                                    bestObjVal = testObjVal;
//                                    logiSeq = neigh_logiSeq;
//                                    improved = true;
//                                }
//                            }
//                            break;
//                        case 1:
//                            if (!feasEndInsertPos.isEmpty()) {
//                                Iterator<int[]> iterator = feasEndInsertPos.iterator();
//                                int[] fea = iterator.next();
//                                feasEndInsertPos.remove(fea);
//                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                                endInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                                if (!checkLambda(neigh_logiSeq)) {
//                                    System.out.println("wrong");
//                                }
//                                double testObjVal = getLogiObj(neigh_logiSeq);
//                                if (testObjVal <= bestObjVal) {
//                                    bestObjVal = testObjVal;
//                                    logiSeq = neigh_logiSeq;
//                                    improved = true;
//                                }
//                            }
//                            break;
//                        case 2:
//                            if (!feasNewInsertPos.isEmpty()) {
//                                Iterator<int[]> iterator = feasNewInsertPos.iterator();
//                                int[] fea = iterator.next();
//                                feasNewInsertPos.remove(fea);
//                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                                newInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                                if (!checkLambda(neigh_logiSeq)) {
//                                    System.out.println("wrong");
//                                }
//                                double testObjVal = getLogiObj(neigh_logiSeq);
//                                if (testObjVal <= bestObjVal) {
//                                    bestObjVal = testObjVal;
//                                    logiSeq = neigh_logiSeq;
//                                    improved = true;
//                                }
//                            }
//                            break;
//                        case 3:
//                            if (!feasSwapPos.isEmpty()) {
//                                Iterator<int[]> iterator = feasSwapPos.iterator();
//                                int[] fea = iterator.next();
//                                feasSwapPos.remove(fea);
//                                ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                                swapOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                                if (!checkLambda(neigh_logiSeq)) {
//                                    System.out.println("wrong");
//                                }
//                                double testObjVal = getLogiObj(neigh_logiSeq);
//                                if (testObjVal <= bestObjVal) {
//                                    bestObjVal = testObjVal;
//                                    logiSeq = neigh_logiSeq;
//                                    improved = true;
//                                }
//                            }
//                            break;
//                    }
//                }
//                iter++;
//            }
//        }
//        if(bestObjVal<targetObj) {
//
//            // 根据物流顺序返回tad的值
//            double[] tadd = new double[logiTaskSet.size()];
//            for (int i = 0; i < logiTaskSet.size(); i++) {
//                tadd[i] = Float.MAX_VALUE;
//            }
//            for (int i = 0; i < logiSeq.size(); i++) {
//                double baseTime = 0;
//                for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                    LogiTask logiTask = para.ItoLogiTask.get(logiSeq.get(i).get(j));
//                    double temp = baseTime + logiTask.getChargingTime() + logiTask.getTransTime() / 2;
//                    if (temp < tadd[logiTask.getId()]) {
//                        tadd[logiTask.getId()] = temp;
//                    }
//                    baseTime = baseTime + logiTask.getProcessTime();
//                }
//            }
//            for (LogiTask logiTask : logiTaskSet) {
//                logiTask.setTad(tadd[logiTask.getId()]);
//            }
////        for (LogiTask logiTask : logiTaskSet) {// 列生成获得的obj基础上减去sum_i (lambda_i*TransTime_(0,i))才是物流子问题的最终目标值
////            logiObj = logiObj - logiTask.getLambda() * (logiTask.getTransTime() / 2);
////        }
////
//            this.logiSeq=listCopy2(logiSeq);
//
//        }
//        return Math.min(bestObjVal,targetObj);
//
//    }

    private boolean checkLambdaProcessTimeRatio1(ArrayList<LogiTask> logiSeq) {
        boolean flag=true;
        if(logiSeq.size()>1){
            for (int j = 0; j < logiSeq.size()-1; j++) {
                double ratio1=logiSeq.get(j).getLambda()/logiSeq.get(j).getProcessTime();
                double ratio2=logiSeq.get(j+1).getLambda()/logiSeq.get(j+1).getProcessTime();
                if(ratio1<ratio2){
                    flag=false;//System.out.println("wrong");
                    System.out.println("veh:"+" ,pos:"+j);
                }
            }
        }
        return flag;
    }
    private boolean checkLambdaProcessTimeRatio(ArrayList<ArrayList<LogiTask>> logiSeq) {
        boolean flag=true;
        for (int i = 0; i < logiSeq.size(); i++) {
            if(logiSeq.get(i).size()>1){
                for (int j = 0; j < logiSeq.get(i).size()-1; j++) {
                    double ratio1=logiSeq.get(i).get(j).getLambda()/logiSeq.get(i).get(j).getProcessTime();
                    double ratio2=logiSeq.get(i).get(j+1).getLambda()/logiSeq.get(i).get(j+1).getProcessTime();
                    if(ratio1<ratio2){
                        flag=false;//System.out.println("wrong");
                        System.out.println("\033[31mveh:\033[0m"+i+" ,pos:"+j);
                    }
                }
            }
        }
        return flag;
    }

    private double solveLogiProblem() throws ClassNotFoundException, IOException, IloException {
//		解物流子问题，更新logiTask的tad
//		System.out.println("\n<<<<<<<<<<< Column Generation For Logistics Problem >>>>>>>>>>>");
        ColGen colGen = new ColGen(para);
        colGen.run();
        double logiObj = colGen.getObj();// 列生成获得的obj是不带 sum_i (lambda_i*TransTime_(0,i))的
//
        optAgvPaths = colGen.rstMasterProblem.getOptAgvPaths();// ArrayList<AgvPath>类型

        // 根据物流顺序返回tad的值
        double[] tad = new double[logiTaskSet.size()];
        for (int i = 0; i < logiTaskSet.size(); i++) {
            tad[i] = Float.MAX_VALUE;
        }
        for (AgvPath agvPath : optAgvPaths) {
            double baseTime = 0;
            for (LogiTask logiTask : agvPath.getExecuteSequence()) {
                double temp = baseTime + logiTask.getChargingTime() + logiTask.getTransTime() / 2;
                if (temp < tad[logiTask.getId()]) {
                    tad[logiTask.getId()] = temp;
                }
//				logiTask.setTad(baseTime +logiTask.getChargingTime()+ logiTask.getTransTime() / 2);
                baseTime = baseTime + logiTask.getProcessTime();
            }
        }
        for (LogiTask logiTask : logiTaskSet) {
            logiTask.setTad(tad[logiTask.getId()]);
        }
        for (LogiTask logiTask : logiTaskSet) {// 列生成获得的obj基础上减去sum_i (lambda_i*TransTime_(0,i))才是物流子问题的最终目标值
            logiObj = logiObj - logiTask.getLambda() * (logiTask.getTransTime() / 2);
        }

        logiSeq.clear();
        for (int r = 0; r < para.numVeh; r++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        for (int r = 0; r < optAgvPaths.size(); r++) {
            AgvPath agvPath=optAgvPaths.get(r);
            for (LogiTask logiTask : agvPath.getExecuteSequence()) {
                logiSeq.get(r).add(logiTask.getId());
            }
        }

        return logiObj;

    }

    public double LagrangHeuristic() {
        String key;
        double timeCur;// 时间游标
        double sum_taf = 0;
        double iniMaxTaf = 0;
        key = prodSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(PMKtoLogiTask.get(key).getTss(), para.moldChangeTime);// 第1个产品第1个机器第1个批次的物料的开始时间（生产问题求得），是整个系统的最早开始时间
        PMKtoLogiTask.get(key).tssHeur_Earl = Math.max(timeCur, PMKtoLogiTask.get(key).getTad());
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).getTad()); // 真正的开始时间是物料到达和最早开工时间取大；
                    PMKtoLogiTask.get(key).tssHeur_Earl = timeCur;
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        para.prodSet.get(prodSeq[p]).tafHeur_Earl = (timeCur);
                    }
                }
            }
//			if (p==para.numProd-1) {
//				iniMaxTaf=timeCur;
//			}
//			timeCur=timeCur+(para.jobshopLength/para.numMach)/para.lineSpeed;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
//		for (int p = 0; p < para.numProd; p++) {
//			tafUB[p]=tafUB[p]+iniMaxTaf-tafUB[para.numProd-1];
//		}
        return sum_taf;
    }

    public double LagrangHeuristic2(int[] prodSeq) {
//		if(Math.random()<0.1) {
//			Random rdm= new Random();
//			int idx1=rdm.nextInt(prodSeq.length);
//			int idx2=rdm.nextInt(prodSeq.length);
//			int temp=prodSeq[idx1];
//			prodSeq[idx1]=prodSeq[idx2];
//			prodSeq[idx2]=temp;
//		}

        logiSeq.clear();
        double[] timeCur = new double[para.numVeh];
        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
//		ArrayList<ArrayList<Integer>> logiSeq= new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq1 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq2 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq3 = new ArrayList<ArrayList<Integer>>();
//		int[] preLogiTaskID=new int[para.numVeh];
        for (int r = 0; r < para.numVeh; r++) {
//			preLogiTaskID[r]=-1;
            logiSeq.add(new ArrayList<Integer>());
            logiSeq1.add(new ArrayList<Integer>());
            logiSeq2.add(new ArrayList<Integer>());
            logiSeq3.add(new ArrayList<Integer>());
        }
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            ArrayList<Integer> idSet = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }

            while (!arrayOfI.isEmpty()) {
                Random random = new Random();
                int index = random.nextInt(arrayOfI.size());
                idSet.add(arrayOfI.get(index));
                arrayOfI.remove(index);
            }

            for (int index = 0; index < idSet.size(); index++) {
//			int machId = arrayofPMKI.get(index)[1];// 对应mach的ID
//			int batchId = arrayofPMKI.get(index)[2];// 对应batch的ID
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                if (Math.random() < prob) {
                    Random random = new Random();
                    min_index = random.nextInt(para.numVeh);

                    logiSeq1.get(min_index).add(logiTaskId);
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                } else {
                    double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    double min_val = Float.MAX_VALUE;
                    for (int r = 0; r < para.numVeh; r++) {
                        if (timeCur[r] + requiredDeliTime < min_val) {
                            min_val = timeCur[r] + requiredDeliTime;
                            min_index = r;
                        }
                    }
                    logiSeq1.get(min_index).add(logiTaskId);
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                }
            }


            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                Random random = new Random();
                min_index = random.nextInt(para.numVeh);

                logiSeq2.get(min_index).add(logiTaskId);
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }


            for (int index = 0; index < idSet.size(); index++) {
//					int machId = arrayofPMKI.get(index)[1];// 对应mach的ID
//					int batchId = arrayofPMKI.get(index)[2];// 对应batch的ID
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
//						Random random2 = new Random();
//						int index2=random2.nextInt(para.numVeh);
//						logiSeq.get(index2).add(logiTaskId);
                int min_index = 0;
                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                double min_val = Float.MAX_VALUE;
                for (int r = 0; r < para.numVeh; r++) {
                    if (timeCur[r] + requiredDeliTime < min_val) {
                        min_val = timeCur[r] + requiredDeliTime;
                        min_index = r;
                    }
                }
                logiSeq3.get(min_index).add(logiTaskId);
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }

        }
        v1 = getObj(prodSeq, logiSeq1);
        v2 = getObj(prodSeq, logiSeq2);
        v3 = getObj(prodSeq, logiSeq3);

        double testObjVal = 0;
        if (v1 <= v2 && v1 <= v3) {
            logiSeq = logiSeq1;
        }
        if (v2 <= v1 && v2 <= v3) {
            logiSeq = logiSeq2;
        }
        if (v3 <= v1 && v3 <= v2) {
            logiSeq = logiSeq3;
        }

        testObjVal = getObj(prodSeq, logiSeq);
        return testObjVal;
    }

    public double LagrangHeuristic3(int[] prodSeq) {
        logiSeq.clear();
        double[] timeCur = new double[para.numVeh];
        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
        int[][] posi1=new int[para.numProd][4];
        int[][] posi2=new int[para.numProd][4];
        int[][] posi3=new int[para.numProd][4];
//		ArrayList<ArrayList<Integer>> logiSeq= new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq1 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq2 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq3 = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> first = new ArrayList<>();
        ArrayList<Integer> last = new ArrayList<>();
//		int[] preLogiTaskID=new int[para.numVeh];
        for (int r = 0; r < para.numVeh; r++) {
//			preLogiTaskID[r]=-1;
            logiSeq.add(new ArrayList<Integer>());
            logiSeq1.add(new ArrayList<Integer>());
            logiSeq2.add(new ArrayList<Integer>());
            logiSeq3.add(new ArrayList<Integer>());
        }
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            ArrayList<Integer> idSet = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }
            first.add(arrayOfI.get(0));
            last.add(arrayOfI.get(arrayOfI.size()-1));

            while (!arrayOfI.isEmpty()) {
                Random random = new Random();
                int index = random.nextInt(arrayOfI.size());
                idSet.add(arrayOfI.get(index));
                arrayOfI.remove(index);
            }

            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                if (Math.random() < prob) {
                    Random random = new Random();
                    min_index = random.nextInt(para.numVeh);

                    logiSeq1.get(min_index).add(logiTaskId);
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                } else {
                    double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    double min_val = Float.MAX_VALUE;
                    for (int r = 0; r < para.numVeh; r++) {
                        if (timeCur[r] + requiredDeliTime < min_val) {
                            min_val = timeCur[r] + requiredDeliTime;
                            min_index = r;
                        }
                    }
                    logiSeq1.get(min_index).add(logiTaskId);
                    if(first.contains(logiTaskId)){
                        posi1[pindex][0]=min_index;
                        posi1[pindex][1]=logiSeq1.get(min_index).size()-1;
                    }
                    if(last.contains(logiTaskId)){
                        posi1[pindex][2]=min_index;
                        posi1[pindex][3]=logiSeq1.get(min_index).size()-1;
                    }
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                }
            }

            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                Random random = new Random();
                min_index = random.nextInt(para.numVeh);

                logiSeq2.get(min_index).add(logiTaskId);
                if(first.contains(logiTaskId)){
                    posi2[pindex][0]=min_index;
                    posi2[pindex][1]=logiSeq2.get(min_index).size()-1;
                }
                if(last.contains(logiTaskId)){
                    posi2[pindex][2]=min_index;
                    posi2[pindex][3]=logiSeq2.get(min_index).size()-1;
                }
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }

            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                double min_val = Float.MAX_VALUE;
                for (int r = 0; r < para.numVeh; r++) {
                    if (timeCur[r] + requiredDeliTime < min_val) {
                        min_val = timeCur[r] + requiredDeliTime;
                        min_index = r;
                    }
                }
                logiSeq3.get(min_index).add(logiTaskId);
                if(first.contains(logiTaskId)){
                    posi3[pindex][0]=min_index;
                    posi3[pindex][1]=logiSeq3.get(min_index).size()-1;
                }
                if(last.contains(logiTaskId)){
                    posi3[pindex][2]=min_index;
                    posi3[pindex][3]=logiSeq3.get(min_index).size()-1;
                }
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }
        }
        v1 = getObj(prodSeq, logiSeq1);
        v2 = getObj(prodSeq, logiSeq2);
        v3 = getObj(prodSeq, logiSeq3);

        int[][] posi=new int[para.numProd][4];
        double testObjVal = 0;
        if (v1 <= v2 && v1 <= v3) {
            logiSeq = logiSeq1;
            posi=posi1;
        }
        if (v2 <= v1 && v2 <= v3) {
            logiSeq = logiSeq2;
            posi=posi2;
        }
        if (v3 <= v1 && v3 <= v2) {
            logiSeq = logiSeq3;
            posi=posi3;
        }

        testObjVal = getObj(prodSeq, logiSeq);
        testObjVal=Math.min(testObjVal,LocalImprove(prodSeq, logiSeq,testObjVal,posi));
        return testObjVal;
    }

    public double LocalSearch_backup(int[] prodSeq) {
        double bestObjVal=Float.MAX_VALUE;

        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
        for (int r = 0; r < para.numVeh; r++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        double[] tad=new double[logiTaskSet.size()];
        double[] tss=new double[logiTaskSet.size()];
        double[] delta=new double[logiTaskSet.size()];
        double[] lambda=new double[logiTaskSet.size()];
        for (LogiTask logiTask:logiTaskSet) {
            lambda[logiTask.getId()]=logiTask.getLambda();
            tad[logiTask.getId()]=logiTask.getTad();
            tss[logiTask.getId()]=logiTask.getTss();
            delta[logiTask.getId()]=tad[logiTask.getId()]-tss[logiTask.getId()];
        }

        ArrayList<ArrayList<Integer>> iList=new ArrayList<ArrayList<Integer>>();
        for (int p = 0; p < prodSeq.length; p++) {
            iList.add(new ArrayList<Integer>());
        }

        for (int p = 0; p < prodSeq.length; p++) {
            ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
            String key;
            double timeCur = 0;
            for (int m = 0; m < para.numMach; m++) {
                if (m > 0) {
                    timeCur = timeCur
                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
                            / para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
                    key = prodSeq[p] + "_" + m + "_" + k;
                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
            for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
                String pmk = tssWithoutLogiConsList.get(s).pmk;
                iList.get(p).add(para.PMKtoLogiTask.get(pmk).getId());
            }
        }

        ArrayList<ArrayList<Integer>> cumb_iList=new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> ini_iList=new ArrayList<ArrayList<Integer>>();
        for (int p = 0; p < prodSeq.length; p++) {
            cumb_iList.add(new ArrayList<Integer>());
            ini_iList.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < iList.size(); i++) {
            for (int j = 0; j < iList.get(i).size(); j++) {
                int idx=iList.get(i).get(j);
                cumb_iList.get(i).add(idx);
                ini_iList.get(i).add(idx);
            }
        }
        int iter=0;
//        while (iter<1000){
//
//            ArrayList<ArrayList<Integer>> neigh_iList= neighbourOperator(cumb_iList,ini_iList);
//            logiSeq=routeGeneration(neigh_iList);
//            double testObjVal = getObj(prodSeq, logiSeq);
//            if(testObjVal<bestObjVal){
//                bestObjVal=testObjVal;
//                for (int i = 0; i < iList.size(); i++) {
//                    for (int j = 0; j < iList.get(i).size(); j++) {
//                        int idx=neigh_iList.get(i).get(j);
//                        cumb_iList.get(i).add(idx);
//                    }
//                }
//            }
//        }

        return bestObjVal;
    }

//    public double LocalSearch(int[] prodSeq) {
//
//        ArrayList<ArrayList<Integer>> logiPerProdList=new ArrayList<ArrayList<Integer>>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            logiPerProdList.add(new ArrayList<Integer>());
//        }
//
//        for (int p = 0; p < prodSeq.length; p++) {
//            ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
//            double timeCur = 0;
//            for (int m = 0; m < para.numMach; m++) {
//                if (m > 0) {
//                    timeCur = timeCur
//                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
//                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
//                            / para.lineSpeed;
//                }
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
//                    String key = prodSeq[p] + "_" + m + "_" + k;
//                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
//                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
//                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//                }
//            }
//            Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
//            for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
//                String pmk = tssWithoutLogiConsList.get(s).pmk;
//                logiPerProdList.get(p).add(para.PMKtoLogiTask.get(pmk).getId());
//            }
//        }
//
//        int[] firstP=new int[para.numProd];
//        int[] lastP=new int[para.numProd];
//        int[] preUP=new int[para.numProd];
//        int[] postUP=new int[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {
//            firstP[p]=logiPerProdList.get(p).get(0);
//            lastP[p]=logiPerProdList.get(p).get( logiPerProdList.get(p).size()-1 );
//            if(p==0){
//                preUP[p]=0;
//                postUP[p]=logiPerProdList.get(p+1).size();
//            }else if(p==para.numProd-1){
//                preUP[p]=logiPerProdList.get(p-1).size();
//                postUP[p]=0;
//            }else {
//                postUP[p]=logiPerProdList.get(p+1).size();
//                preUP[p]=logiPerProdList.get(p-1).size();
//            }
//        }
//
//
//        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
//        for (int r = 0; r < para.numVeh; r++) {
//            logiSeq.add(new ArrayList<Integer>());
//        }
//        double[] tad=new double[logiTaskSet.size()];
//        double[] tss=new double[logiTaskSet.size()];
//        double[] delta=new double[logiTaskSet.size()];
//        double[] lambda=new double[logiTaskSet.size()];
//        for (LogiTask logiTask:logiTaskSet) {
//            lambda[logiTask.getId()]=logiTask.getLambda();
//            tad[logiTask.getId()]=logiTask.getTad();
//            tss[logiTask.getId()]=logiTask.getTss();
//            delta[logiTask.getId()]=tad[logiTask.getId()]-tss[logiTask.getId()];
//        }
//
////        ArrayList<Integer> iList=new ArrayList<Integer>();
////        ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
////        double timeCur = 0;
////        String key;
////        for (int p = 0; p < prodSeq.length; p++) {
////            for (int m = 0; m < para.numMach; m++) {
////                if (m > 0) {
////                    timeCur = timeCur
////                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
////                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
////                            / para.lineSpeed;
////                }
////                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
////                    key = prodSeq[p] + "_" + m + "_" + k;
////                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
////                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
////                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
////                }
////            }
////            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
////                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
////                    / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
////        }
////        Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
////        for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
////            String pmk = tssWithoutLogiConsList.get(s).pmk;
////            iList.add(para.PMKtoLogiTask.get(pmk).getId());
////        }
//
////        int[] iList2=new int[logiTaskSet.size()];
////        tssWithoutLogiConsList.clear();
////        for (LogiTask logiTask:logiTaskSet) {
////            key=logiTask.getProd().getId() + "_" + logiTask.getMach().getId() + "_" + logiTask.getBatchId();
////            tssWithoutLogiConsList.add(new TssWithoutLogiCons(key, logiTask.getLambda()));
////        }
////        Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
////        for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
////            String pmk = tssWithoutLogiConsList.get(s).pmk;
////            iList2[s]=(para.PMKtoLogiTask.get(pmk).getId());
////        }
//
//        ArrayList<ArrayList<Integer>> idSetPerProd=new ArrayList<ArrayList<Integer>>();
//        for (int p = 0; p < para.numProd; p++) {
//            idSetPerProd.add(new ArrayList<Integer>());
//        }
//        for (LogiTask logiTask: para.logiTaskSet) {
//            idSetPerProd.get(logiTask.getProd().getId()).add(logiTask.getId());
//        }
//        ArrayList<Integer> iList=new ArrayList<Integer>();
//        for (int p = 0; p < para.numProd; p++) {
//            for (int i = 0; i < idSetPerProd.get(prodSeq[p]).size(); i++) {
//                iList.add(idSetPerProd.get(prodSeq[p]).get(i));
//            }
//        }
//
//        ArrayList<Integer> cumb_iList=listCopy(iList);
////        ArrayList<Integer> ini_iList=listCopy(iList);
//        double bestObjVal=Float.MAX_VALUE;
//        int iter=0;
//        while (iter<3000){
//            ArrayList<Integer> neigh_iList=listCopy(cumb_iList);
//            ArrayList<ArrayList<Integer>> idxPerProd=new ArrayList<ArrayList<Integer>>();
//            for (int i = 0; i < para.numProd; i++) {
//                idxPerProd.add(new ArrayList<Integer>());
//            }
//            for (int i = 0; i < neigh_iList.size(); i++) {
//                int p=para.ItoLogiTask.get(neigh_iList.get(i)).getProd().getId();
//                idxPerProd.get(p).add(i);
//            }
//            Random rnd = new Random();
//            int origPosition=rnd.nextInt(neigh_iList.size());
//            int rndLogiTaskId=neigh_iList.get(origPosition);
//            int rndP=para.ItoLogiTask.get(rndLogiTaskId).getProd().getId();
//            int posiLB,posiUB,newPosition;
////            if (Arrays.binarySearch(firstP,rndLogiTaskId)>=0){
////                posiLB=Math.max(origPosition-preUP[rndP],0);
////                ArrayList<Integer> candi=new ArrayList<Integer>();
////                for (int i = posiLB; i <origPosition ; i++) {
////                    candi.add(neigh_iList.get(i));
////                }
////                for (int i = 0; i < idxPerProd.get(rndP).size(); i++) {
////                    if (!candi.contains(idxPerProd.get(rndP).get(i))){
////                        candi.add(idxPerProd.get(rndP).get(i));
////                    }
////                }
////                newPosition=candi.get(rnd.nextInt(candi.size()));//rnd.nextInt(posiUB+1-posiLB)+posiLB;
////            } else if (Arrays.binarySearch(lastP,rndLogiTaskId)>=0) {
////                posiUB=Math.min(origPosition+postUP[rndP],neigh_iList.size()-1);
////                ArrayList<Integer> candi=new ArrayList<Integer>();
////                for (int i = origPosition; i <posiUB+1 ; i++) {
////                    candi.add(neigh_iList.get(i));
////                }
////                for (int i = 0; i < idxPerProd.get(rndP).size(); i++) {
////                    if (!candi.contains(idxPerProd.get(rndP).get(i))){
////                        candi.add(idxPerProd.get(rndP).get(i));
////                    }
////                }
////                newPosition=candi.get(rnd.nextInt(candi.size()));
////            }else {
////                newPosition=idxPerProd.get(rndP).get(rnd.nextInt(idxPerProd.get(rndP).size()));//rnd.nextInt(posiUB+1-posiLB)+posiLB;
////            }
////            Collections.swap(neigh_iList,origPosition,newPosition);
//            logiSeq=routeGeneration0(neigh_iList);
//            double testObjVal = getObj(prodSeq, logiSeq);
//            if(testObjVal<=bestObjVal){
//                bestObjVal=testObjVal;
//                cumb_iList=listCopy1(neigh_iList);
//            }
//            iter++;
//        }
//
//        return bestObjVal;
//    }

    public double LocalImprove(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq, double obj, int[][] posi) {
        //posi[p][a]记录了产品p的第一个和最后一个logiTaskId所在的车辆a1和位置a2
        double new_obj=obj;

        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
        for (int p = 0; p < prodSeq.length; p++) {
            prodPriority.add(prodSeq[p]);
        }
        for (int p = 0; p < prodSeq.length; p++) {
            int ovfirst = posi[p][0];
            int opfirst = posi[p][1];
            int ovlast = posi[p][2];
            int oplast = posi[p][3];
            if(prodSeq[0] == p){//产品p是第一个
                HashSet<int[]> insertSet = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                int o_LogiTaskId = logiSeq.get(ovlast).get(oplast);
                int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
                int o_Priority = prodPriority.indexOf(o_ProdId);
                for (int i = 0; i < logiSeq.size(); i++) {
                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
                        int[] fea = {i, j};
                        int n_LogiTaskId = logiSeq.get(i).get(j);
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o_Priority+1) {
                            insertSet.add(fea);
                        }
                    }
                }
                Iterator<int[]> it = insertSet.iterator();
                while(it.hasNext()){
                    int[] fea=it.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovlast,oplast,fea[0],fea[1]);
//                    it.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                    }
                }
            } else if (prodSeq[prodSeq.length-1] == p) {//产品p是最后一个
                HashSet<int[]> insertSet = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                int o_LogiTaskId = logiSeq.get(ovfirst).get(opfirst);
                int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
                int o_Priority = prodPriority.indexOf(o_ProdId);
                for (int i = 0; i < logiSeq.size(); i++) {
                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
                        int[] fea = {i, j};
                        int n_LogiTaskId = logiSeq.get(i).get(j);
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o_Priority-1) {
                            insertSet.add(fea);
                        }
                    }
                }
                Iterator<int[]> it = insertSet.iterator();
                while(it.hasNext()){
                    int[] fea=it.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovfirst,opfirst,fea[0],fea[1]);
//                    it.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                    }
                }
            } else {
                HashSet<int[]> insertSet1 = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                HashSet<int[]> insertSet2 = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                int o1_LogiTaskId = logiSeq.get(ovfirst).get(opfirst);
                int o1_ProdId = para.ItoLogiTask.get(o1_LogiTaskId).getProd().getId();
                int o1_Priority = prodPriority.indexOf(o1_ProdId);
                int o2_LogiTaskId = logiSeq.get(ovlast).get(oplast);
                int o2_ProdId = para.ItoLogiTask.get(o2_LogiTaskId).getProd().getId();
                int o2_Priority = prodPriority.indexOf(o2_ProdId);
                for (int i = 0; i < logiSeq.size(); i++) {
                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
                        int[] fea = {i, j};
                        int n_LogiTaskId = logiSeq.get(i).get(j);
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o1_Priority-1) {
                            insertSet1.add(fea);
                        }
                        if (n_Priority == o2_Priority+1) {
                            insertSet2.add(fea);
                        }
                    }
                }
                Iterator<int[]> it1 = insertSet1.iterator();
                while(it1.hasNext()){
                    int[] fea=it1.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovfirst,opfirst,fea[0],fea[1]);
//                    it1.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                    }
                }
                Iterator<int[]> it2 = insertSet2.iterator();
                while(it2.hasNext()){
                    int[] fea=it2.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovlast,oplast,fea[0],fea[1]);
//                    it2.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                    }
                }
            }
        }

        return new_obj;
    }
//    public double LocalSearch(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq) {
//
////        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
////        for (int p = 0; p < para.numVeh; p++) {
////            logiSeq.add(new ArrayList<Integer>());
////        }
////
////        double[] time=new double[para.numVeh];
////        for (int p = 0; p < prodSeq.length; p++) {
////            ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
////            String key;
////            double timeCur = 0;
////            for (int m = 0; m < para.numMach; m++) {
////                if (m > 0) {
////                    timeCur = timeCur
////                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
////                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
////                            / para.lineSpeed;
////                }
////                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
////                    key = prodSeq[p] + "_" + m + "_" + k;
////                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
////                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
////                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
////                }
////            }
////            Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
////
////
////            for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
////                String pmk = tssWithoutLogiConsList.get(s).pmk;
//////                Random rnd=new Random();
//////                logiSeq.get(rnd.nextInt(para.numVeh)).add(para.PMKtoLogiTask.get(pmk).getId());
////                int logiTaskId=para.PMKtoLogiTask.get(pmk).getId();
////                int min_index = 0;
////                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
////                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
////                double min_val = Float.MAX_VALUE;
////                for (int r = 0; r < para.numVeh; r++) {
////                    if (time[r] + requiredDeliTime < min_val) {
////                        min_val = time[r] + requiredDeliTime;
////                        min_index = r;
////                    }
////                }
////                logiSeq.get(min_index).add(logiTaskId);
////                time[min_index] = time[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
////            }
////        }
//
//
//        double bestObjVal=getObj(prodSeq, logiSeq);
//        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            prodPriority.add(prodSeq[p]);
//        }
//
//        ArrayList<ArrayList<Integer>> best_logiSeq = listCopy2(logiSeq);
//
//        int iter=0;
//        Random rnd=new Random();
//        while (iter<1000) {
//            //随机找一个位置
//            //根据该位置返回可行的替换位置
//            //将两个位置做替换
//
//            int ov = rnd.nextInt(para.numVeh);
//            int op = rnd.nextInt(best_logiSeq.get(ov).size());
//            int o_LogiTaskId = best_logiSeq.get(ov).get(op);
//            int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
//            int o_Priority = prodPriority.indexOf(o_ProdId);
//
//            HashSet<int[]> feasInnerInsertPos = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
//            HashSet<int[]> feasOuterInsertPos = new HashSet<>();//允许跨优先级的插入（插入在指定位置之前）
//            HashSet<int[]> feasAfterInsertPos = new HashSet<>();//插入在指定位置之后
//            HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
//            HashSet<int[]> feasSwapPos = new HashSet<>();//相同优先级的交换
//
//            for (int i = 0; i < best_logiSeq.size(); i++) {
//                for (int j = 0; j < best_logiSeq.get(i).size(); j++) {
//                    int[] fea = {i, j};
//                    int n_LogiTaskId = best_logiSeq.get(i).get(j);
//                    if (o_LogiTaskId != n_LogiTaskId) {
//                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
//                        int n_Priority = prodPriority.indexOf(n_ProdId);
//                        if (n_Priority == o_Priority) {
//                            feasSwapPos.add(fea);
//                            feasInnerInsertPos.add(fea);
//                        }
//                        if (n_Priority < o_Priority) {
//                            feasOuterInsertPos.add(fea);
//                        }
//                        if (j < best_logiSeq.get(i).size() - 1) {
//                            int n2_LogiTaskId = best_logiSeq.get(i).get(j + 1);
//                            int n2_ProdId = para.ItoLogiTask.get(n2_LogiTaskId).getProd().getId();
//                            int n2_Priority = prodPriority.indexOf(n2_ProdId);
//                            if (n_Priority <= o_Priority && n2_Priority > o_Priority) {
//                                feasAfterInsertPos.add(fea);
//                            }
//                        } else {
//                            if (n_Priority <= o_Priority) {
//                                feasEndInsertPos.add(fea);
//                            }
//                        }
//                    }
//                }
//            }
//
//            //找一个最好的领域
//            ArrayList<ArrayList<Integer>> bestNeigh_logiSeq=new ArrayList<ArrayList<Integer>>();
//            Random rndOperator = new Random();
//            boolean improve=false;
//            while(!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty() || !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty()){
//                switch (rndOperator.nextInt(4)) {
//                    case 0:
//                        if (!feasInnerInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasInnerInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasInnerInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            innerInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                    case 1:
//                        if (!feasAfterInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasAfterInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasAfterInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            afterInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                    case 2:
//                        if (!feasEndInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasEndInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasEndInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            endInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                    case 3:
//                        if (!feasSwapPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasSwapPos.iterator();
//                            int[] fea = iterator.next();
//                            feasSwapPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            swapOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                }
//
//            }
//            for (int[] fea : feasOuterInsertPos) {
//                ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                outerInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                if (testObjVal <= bestObjVal) {
//                    bestObjVal = testObjVal;
//                    bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                    improve=true;
//                }
//            }
//
//
////            swapOperator(neigh_logiSeq,ov,op,nv,np);
////            afterInsertOperator(neigh_logiSeq,ov,op,nv,np);
////            endInsertOperator(neigh_logiSeq,ov,op,nv,np);
////            innerInsertOperator(neigh_logiSeq,ov,op,nv,np);
////            outerInsertOperator(neigh_logiSeq,ov,op,nv,np);
//            if(improve) {
//                best_logiSeq = bestNeigh_logiSeq;
//            }
//
//            iter++;
//        }
//
//        return bestObjVal;
//    }

    public double LocalSearch(int[] prodSeq,  ArrayList<ArrayList<Integer>> ini_logiSeq) {

        if(Math.random()<0.1){
            Random rnd=new Random();
            int p_change=rnd.nextInt(para.numProd-1);
            int temp = prodSeq[p_change+1];
            prodSeq[p_change+1]=prodSeq[p_change];
            prodSeq[p_change]=temp;
        }
        ArrayList<ArrayList<Integer>> cumb_logiSeq=new ArrayList<ArrayList<Integer>>();
        for (int r = 0; r < para.numVeh; r++) {
            cumb_logiSeq.add(new ArrayList<Integer>());
        }

        ArrayList<Integer> first = new ArrayList<>();
        ArrayList<Integer> last = new ArrayList<>();
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }
            first.add(arrayOfI.get(0));
            last.add(arrayOfI.get(arrayOfI.size() - 1));
        }

        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
        for (int p = 0; p < prodSeq.length; p++) {
            prodPriority.add(prodSeq[p]);
        }
        int[] priorityLogiTask=new int[logiTaskSet.size()];
        for (int i = 0; i < ini_logiSeq.size(); i++) {
            for (int j = 0; j < ini_logiSeq.get(i).size(); j++) {
                int logiTaskId=ini_logiSeq.get(i).get(j);
                int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
                priorityLogiTask[logiTaskId]=prodPriority.indexOf(prodId);
            }
        }
        for (int i = 0; i < ini_logiSeq.size(); i++) {
            ArrayList<LogiTaskPrior> priorSet=new ArrayList<LogiTaskPrior>();
            for (int j = 0; j < ini_logiSeq.get(i).size(); j++) {
                int logiTaskId = ini_logiSeq.get(i).get( j );
                priorSet.add(new LogiTaskPrior(logiTaskId,priorityLogiTask[logiTaskId]));
//                System.out.print(priorityLogiTask[logiTaskId]+" ");
            }
//            System.out.println(priorSet.size());
            Collections.sort(priorSet, new LogiTaskSortByPrior());

            for (int j = 0; j < priorSet.size(); j++) {
                cumb_logiSeq.get(i).add(priorSet.get(j).id);
            }

        }

        double bestObjVal=Math.min(getObj(prodSeq, cumb_logiSeq),getObj(prodSeq, ini_logiSeq));

        ArrayList<ArrayList<Integer>> logiSeq = listCopy2(cumb_logiSeq);

        int iter=0;
        Random rnd=new Random();
        while (iter<2000) {
            int ov;
            int op ;
            while (true){
                int a = rnd.nextInt(para.numVeh);
                if (logiSeq.get(a).size()>1){
                    ov=a;
                    op=rnd.nextInt(logiSeq.get(a).size());
                    break;
                }
            }
            int o_LogiTaskId = logiSeq.get(ov).get(op);
            int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
            int o_Priority = prodPriority.indexOf(o_ProdId);

            HashSet<int[]> feasInnerInsertPos = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
            HashSet<int[]> feasAfterInsertPos = new HashSet<>();//插入在指定位置之后
            HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
            HashSet<int[]> feasSwapPos = new HashSet<>();//相同优先级的交换

            for (int i = 0; i < logiSeq.size(); i++) {
                for (int j = 0; j < logiSeq.get(i).size(); j++) {
                    int[] fea = {i, j};
                    int n_LogiTaskId = logiSeq.get(i).get(j);
                    if (o_LogiTaskId != n_LogiTaskId) {
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o_Priority) {
                            feasSwapPos.add(fea);
                            feasInnerInsertPos.add(fea);
                        }
                        if (j < logiSeq.get(i).size() - 1) {
                            int n2_LogiTaskId = logiSeq.get(i).get(j + 1);
                            int n2_ProdId = para.ItoLogiTask.get(n2_LogiTaskId).getProd().getId();
                            int n2_Priority = prodPriority.indexOf(n2_ProdId);
                            if (n_Priority <= o_Priority && n2_Priority > o_Priority) {
                                feasAfterInsertPos.add(fea);
                            }
                        } else {
                            if (n_Priority <= o_Priority) {
                                feasEndInsertPos.add(fea);
                            }
                        }
                    }
                }
            }
            Random rndOperator = new Random();
            boolean improved=false;
            while( !improved && (!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty() || !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty() ) ){
                switch (rndOperator.nextInt(4)) {
                    case 0:
                        if (!feasInnerInsertPos.isEmpty()) {
                            Iterator<int[]> iterator = feasInnerInsertPos.iterator();
                            int[] fea = iterator.next();
                            feasInnerInsertPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            innerInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                logiSeq= neigh_logiSeq;
//                                if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                    int a=0;
//                                }
                                improved = true;
                            }
                        }
                        break;
                    case 1:
                        if (!feasEndInsertPos.isEmpty()) {
                            Iterator<int[]> iterator = feasEndInsertPos.iterator();
                            int[] fea = iterator.next();
                            feasEndInsertPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            endInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                logiSeq= neigh_logiSeq;
                                improved = true;
                            }
                        }
                        break;
                    case 2:
                        if (!feasAfterInsertPos.isEmpty()) {
                            Iterator<int[]> iterator = feasAfterInsertPos.iterator();
                            int[] fea = iterator.next();
                            feasAfterInsertPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            afterInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                logiSeq= neigh_logiSeq;
                                improved = true;
                            }
                        }
                        break;
                    case 3:
                        if (!feasSwapPos.isEmpty()) {
                            Iterator<int[]> iterator = feasSwapPos.iterator();
                            int[] fea = iterator.next();
                            feasSwapPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            swapOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                logiSeq= neigh_logiSeq;
                                improved = true;
                            }
                        }
                        break;
                }
                //做完简单的领域搜索后再试探一下跨优先级的更优解，但是这个解不会被记录下来
//                int[][] posi=new int[para.numProd][4];
//                for (int i = 0; i < logiSeq.size(); i++) {
//                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                        int logiTaskId=logiSeq.get(i).get(j);
//                        int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
//                        if(first.contains(logiTaskId)){
//                            posi[prodId][0]=i;
//                            posi[prodId][1]=j;
//                        }
//                        if(last.contains(logiTaskId)){
//                            posi[prodId][2]=i;
//                            posi[prodId][3]=j;
//                        }
//                    }
//                }
//                ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
//                double testObjVal=LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi);
//                if (testObjVal <= bestObjVal) {
//                    bestObjVal = testObjVal;
//                    improved = true;
//                }
            }
            iter++;
        }

        int[][] posi=new int[para.numProd][4];
        for (int i = 0; i < logiSeq.size(); i++) {
            for (int j = 0; j < logiSeq.get(i).size(); j++) {
                int logiTaskId=logiSeq.get(i).get(j);
                int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
                if(first.contains(logiTaskId)){
                    posi[prodId][0]=i;
                    posi[prodId][1]=j;
                }
                if(last.contains(logiTaskId)){
                    posi[prodId][2]=i;
                    posi[prodId][3]=j;
                }
            }
        }
        ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
        bestObjVal=Math.min(bestObjVal,LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi));

        return bestObjVal;
    }

//    public double LocalSearch(int[] prodSeq) {
//
//
//        ArrayList<Integer> first = new ArrayList<>();
//        ArrayList<Integer> last = new ArrayList<>();
//        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
//            int pindex = prodSeq[p];
//            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
//            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
//            for (int index = 0; index < arrayofPMKI.size(); index++) {
//                arrayOfI.add(arrayofPMKI.get(index)[3]);
//            }
//            first.add(arrayOfI.get(0));
//            last.add(arrayOfI.get(arrayOfI.size() - 1));
//        }
//
//        ArrayList<ArrayList<Integer>> ini_logiSeq= new ArrayList<ArrayList<Integer>>();
//        for (int i = 0; i < para.numVeh; i++) {
//            ini_logiSeq.add(new ArrayList<Integer>());
//        }
//        for (int p = 0; p < prodSeq.length; p++) {
//            int pindex = prodSeq[p];
//            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
//            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
//            ArrayList<Integer> idSet = new ArrayList<Integer>();
//            for (int index = 0; index < arrayofPMKI.size(); index++) {
//                arrayOfI.add(arrayofPMKI.get(index)[3]);
//            }
//            for (int s = 0; s < arrayOfI.size(); s++) {
//                Random rnd=new Random();
//                ini_logiSeq.get(rnd.nextInt(para.numVeh)).add(arrayOfI.get(s));
//            }
//        }
//
//        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            prodPriority.add(prodSeq[p]);
//        }
//
//        double bestObjVal=getObj(prodSeq, ini_logiSeq);
//
//        ArrayList<ArrayList<Integer>> logiSeq = listCopy2(ini_logiSeq);
//
//        int iter=0;
//        Random rnd=new Random();
//        while (iter<1000) {
//            int ov;
//            int op ;
//            while (true){
//                int a = rnd.nextInt(para.numVeh);
//                if (logiSeq.get(a).size()>1){
//                    ov=a;
//                    op=rnd.nextInt(logiSeq.get(a).size());
//                    break;
//                }
//            }
//            int o_LogiTaskId = logiSeq.get(ov).get(op);
//            int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
//            int o_Priority = prodPriority.indexOf(o_ProdId);
//
//            HashSet<int[]> feasInnerInsertPos = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
//            HashSet<int[]> feasAfterInsertPos = new HashSet<>();//插入在指定位置之后
//            HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
//            HashSet<int[]> feasSwapPos = new HashSet<>();//相同优先级的交换
//
//            for (int i = 0; i < logiSeq.size(); i++) {
//                for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                    int[] fea = {i, j};
//                    int n_LogiTaskId = logiSeq.get(i).get(j);
//                    if (o_LogiTaskId != n_LogiTaskId) {
//                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
//                        int n_Priority = prodPriority.indexOf(n_ProdId);
//                        if (n_Priority == o_Priority) {
//                            feasSwapPos.add(fea);
//                            feasInnerInsertPos.add(fea);
//                        }
//                        if (j < logiSeq.get(i).size() - 1) {
//                            int n2_LogiTaskId = logiSeq.get(i).get(j + 1);
//                            int n2_ProdId = para.ItoLogiTask.get(n2_LogiTaskId).getProd().getId();
//                            int n2_Priority = prodPriority.indexOf(n2_ProdId);
//                            if (n_Priority <= o_Priority && n2_Priority > o_Priority) {
//                                feasAfterInsertPos.add(fea);
//                            }
//                        } else {
//                            if (n_Priority <= o_Priority) {
//                                feasEndInsertPos.add(fea);
//                            }
//                        }
//                    }
//                }
//            }
//            Random rndOperator = new Random();
//            boolean improved=false;
//            while( !improved && (!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty() || !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty() ) ){
//                switch (rndOperator.nextInt(4)) {
//                    case 0:
//                        if (!feasInnerInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasInnerInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasInnerInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            innerInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                    int a=0;
//                                }
//                                improved = true;
//                            }
//                        }
//                        break;
//                    case 1:
//                        if (!feasEndInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasEndInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasEndInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            endInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                improved = true;
//                            }
//                        }
//                        break;
//                    case 2:
//                        if (!feasAfterInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasAfterInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasAfterInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            afterInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                improved = true;
//                            }
//                        }
//                        break;
//                    case 3:
//                        if (!feasSwapPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasSwapPos.iterator();
//                            int[] fea = iterator.next();
//                            feasSwapPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            swapOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                improved = true;
//                            }
//                        }
//                        break;
//                }
//                //做完简单的领域搜索后再试探一下跨优先级的更优解，但是这个解不会被记录下来
//                int[][] posi=new int[para.numProd][4];
//                for (int i = 0; i < logiSeq.size(); i++) {
//                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                        int logiTaskId=logiSeq.get(i).get(j);
//                        int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
//                        if(first.contains(logiTaskId)){
//                            posi[prodId][0]=i;
//                            posi[prodId][1]=j;
//                        }
//                        if(last.contains(logiTaskId)){
//                            posi[prodId][2]=i;
//                            posi[prodId][3]=j;
//                        }
//                    }
//                }
//                ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
//                double testObjVal=LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi);
//                if (testObjVal <= bestObjVal) {
//                    bestObjVal = testObjVal;
//                    improved = true;
//                }
//            }
//            iter++;
//        }
//
//        int[][] posi=new int[para.numProd][4];
//        for (int i = 0; i < logiSeq.size(); i++) {
//            for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                int logiTaskId=logiSeq.get(i).get(j);
//                int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
//                if(first.contains(logiTaskId)){
//                    posi[prodId][0]=i;
//                    posi[prodId][1]=j;
//                }
//                if(last.contains(logiTaskId)){
//                    posi[prodId][2]=i;
//                    posi[prodId][3]=j;
//                }
//            }
//        }
//        ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
//        bestObjVal=Math.min(bestObjVal,LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi));
//
//        return bestObjVal;
//    }

//    public boolean checkPriority(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq){
//        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            prodPriority.add(prodSeq[p]);
//        }
//        boolean flag=true;
//        for (int i = 0; i < logiSeq.size(); i++) {
//            if (logiSeq.get(i).size()>1) {
//                for (int j = 0; j < logiSeq.get(i).size() - 1; j++) {
//                    int prior1=prodPriority.indexOf(  para.ItoLogiTask.get(  logiSeq.get(i).get(j)  ).getProd().getId() );
//                    int prior2=prodPriority.indexOf(  para.ItoLogiTask.get(  logiSeq.get(i).get(j+1)  ).getProd().getId() );
//                    if (prior1>prior2){
//                        flag=false;
//                        System.out.println(i+":  "+j+ " and " +(j+1));
//                    }
//                }
//            }
//        }
//        return flag;
//    }

    public void insertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        if(original_position<=new_position) {
            logiSeq.get(new_veh).add(new_position, temp1);
            logiSeq.get(original_veh).remove(original_position);
        }else{
            logiSeq.get(original_veh).remove(original_position);
            logiSeq.get(new_veh).add(new_position, temp1);
        }
    }
    public void swapOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        int temp2=logiSeq.get(new_veh).get(new_position);
        logiSeq.get(original_veh).set(original_position,temp2);
        logiSeq.get(new_veh).set(new_position,temp1);
    }
    public void innerInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        if(original_position<=new_position) {
            logiSeq.get(new_veh).add(new_position, temp1);
            logiSeq.get(original_veh).remove(original_position);
        }else{
            logiSeq.get(original_veh).remove(original_position);
            logiSeq.get(new_veh).add(new_position, temp1);
        }
    }
    public void beforeInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        logiSeq.get(original_veh).remove(original_position);
        logiSeq.get(new_veh).add(new_position, temp1);
    }
    public void newInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        logiSeq.get(original_veh).remove(original_position);
        logiSeq.get(new_veh).add(temp1);
    }
//    public void outerInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
//        int temp1=logiSeq.get(original_veh).get(original_position);
//        logiSeq.get(original_veh).remove(original_position);
//        logiSeq.get(new_veh).add(new_position,temp1);
//    }
    public void afterInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        logiSeq.get(new_veh).add(new_position+1,temp1);
        logiSeq.get(original_veh).remove(original_position);
    }
    public void endInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        logiSeq.get(new_veh).add(temp1);
        logiSeq.get(original_veh).remove(original_position);
    }


//    public ArrayList<Integer> listCopy1(ArrayList<Integer> a){
//        ArrayList<Integer> b=new ArrayList<Integer>();
//        for (int i = 0; i < a.size(); i++) {
//            b.add(a.get(i));
//        }
//        return b;
//    }
    public ArrayList<ArrayList<Integer>> listCopy2(ArrayList<ArrayList<Integer>> a){
        ArrayList<ArrayList<Integer>> b=new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < a.size(); i++) {
            b.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < a.get(i).size(); j++) {
                b.get(i).add(a.get(i).get(j));
            }
        }
        return b;
    }

//    public ArrayList<ArrayList<Integer>> swapOperator(ArrayList<ArrayList<Integer>> cumb_iList){
//        //initialize neigh_iList
//        ArrayList<ArrayList<Integer>> neigh_iList=new ArrayList<ArrayList<Integer>>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            neigh_iList.add(new ArrayList<Integer>());
//        }
//        for (int i = 0; i < cumb_iList.size(); i++) {
//            for (int j = 0; j < cumb_iList.get(i).size(); j++) {
//                int idx=cumb_iList.get(i).get(j);
//                neigh_iList.get(i).add(idx);
//            }
//        }
//
//
//        Random rnd=new Random();
//        int rndP=rnd.nextInt(para.numProd);
//        int rnd1=rnd.nextInt(neigh_iList.get(rndP).size());
//        int rnd2=rnd.nextInt(neigh_iList.get(rndP).size());
//        Collections.swap(neigh_iList.get(rndP),rnd1,rnd2);
//
//        return neigh_iList;
//    }

//    public ArrayList<ArrayList<Integer>> listCopy(ArrayList<ArrayList<Integer>> list){
//
//    }

//    public ArrayList<ArrayList<Integer>> routeGeneration0(ArrayList<Integer> iList){
//        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
//        for (int r = 0; r < para.numVeh; r++) {
//            logiSeq.add(new ArrayList<Integer>());
//        }
//        double[] timeCur = new double[para.numVeh];
//        for (int index = 0; index < iList.size(); index++) {
//            int logiTaskId = iList.get(index);// 对应logiTask的ID
//            int min_index = 0;
//            double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
//                    + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
//            double min_val = Float.MAX_VALUE;
//            for (int r = 0; r < para.numVeh; r++) {
//                if (timeCur[r] + requiredDeliTime < min_val) {
//                    min_val = timeCur[r] + requiredDeliTime;
//                    min_index = r;
//                }
//            }
//            logiSeq.get(min_index).add(logiTaskId);
//            timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
//        }
//        return logiSeq;
//    }
//    public ArrayList<ArrayList<Integer>> routeGeneration(ArrayList<ArrayList<Integer>> iList){
//        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
//        for (int r = 0; r < para.numVeh; r++) {
//            logiSeq.add(new ArrayList<Integer>());
//        }
//        double[] timeCur = new double[para.numVeh];
//        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
//            ArrayList<Integer> idSet=iList.get(p);
//
//            for (int index = 0; index < idSet.size(); index++) {
//                int logiTaskId = idSet.get(index);// 对应logiTask的ID
//                int min_index = 0;
//                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
//                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
//                double min_val = Float.MAX_VALUE;
//                for (int r = 0; r < para.numVeh; r++) {
//                    if (timeCur[r] + requiredDeliTime < min_val) {
//                        min_val = timeCur[r] + requiredDeliTime;
//                        min_index = r;
//                    }
//                }
//                logiSeq.get(min_index).add(logiTaskId);
//                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
//            }
//        }
//        return logiSeq;
//    }
//    public void getTad(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq) {
//        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
//            double cur = 0;
//            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
//                para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur = cur
//                        + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
//                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
//                cur = para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur;
//                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
//                        / para.vehSpeed;
//            }
//        }
//
//        String key;
//        double timeCur = 0;// 时间游标
//        double sum_taf = 0;
//        double[] taf = new double[para.numProd];
//        key = prodSeq[0] + "_" + 0 + "_" + 0;
//        timeCur = PMKtoLogiTask.get(key).tadHeur;// 第1个产品第1个机器第1个批次的物料的开始时间 必然是其对应的tadHeur(已经算好了)
//        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1);
//                }
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    key = prodSeq[p] + "_" + m + "_" + k;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
//                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//                }
//            }
//            taf[p] = timeCur;
//            sum_taf = sum_taf + timeCur;
//        }
//
//        for (int p = 0; p < para.numProd; p++) {
//            System.out.print("-> " + prodSeq[p] + "(" + taf[p] + ") ");
//        }
//        System.out.println();
//
//        for (int i = 0; i < logiSeq.size(); i++) {
//            for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                System.out.print("-> " + logiSeq.get(i).get(j) + "("
//                        + para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur + ") ");
//            }
//            System.out.println();
//        }
//        System.out.println();
//
//    }

    public void updateSubgradient() {
//		update subgradient with the solution of lagrangian relaxation problem given the lagrangian multipliers
        for (LogiTask logiTask : logiTaskSet) {
//			logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()+earliestReadyTime));
            logiTask.setSubgrad(logiTask.getTad() - logiTask.getTss());
//			logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()-makespan));
//			logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()+para.getDistance(para.depot, para.machSet.get(0))*2*para.consumingRate*para.chargingRate+para.getDistance(para.depot, para.machSet.get(0))/para.vehSpeed ));
        }

//		int cnt=0;
//		for(LogiTask logiTask : logiTaskSet) {
//			if (logiTask.getSubgrad()<=0) {
//				cnt++;
//			}
//		}
//		if (cnt==logiTaskSet.size()) {
//			for (LogiTask logiTask : logiTaskSet) {
//				logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()));
//			}
//		}
    }

    public void updateStepsize(double surDual) {
        if (improved) {
            iter_for_stepsize_factor = 0;
        } else {
            iter_for_stepsize_factor++;
        }
        if (iter_for_stepsize_factor >= max_iter_no_improve_for_stepsize_factor) {// 若目标值在固定次数的迭代后仍无改进，stepsizeFactor值将减半
            stepsizeFactor = stepsizeFactor / 2;
            iter_for_stepsize_factor = 0;
        }
        double numerator = stepsizeFactor * (bestObjUpBound - surDual);
        double denominator = 0;
        for (LogiTask logiTask : logiTaskSet) {
            denominator = denominator
                    + ((logiTask.getTad() - logiTask.getTss()) * (logiTask.getTad() - logiTask.getTss()));
        }
        stepsize = numerator / denominator;
    }

    public void updateStepsize() {
        boolean bb = true;
        if (bb) {
//			stepsizeFactor = 0.02;
            if (improved) {
                iter_for_stepsize_factor = 0;
            } else {
                iter_for_stepsize_factor++;
            }
            if (iter_for_stepsize_factor >= max_iter_no_improve_for_stepsize_factor) {// 若目标值在固定次数的迭代后仍无改进，stepsizeFactor值将减半
                stepsizeFactor = stepsizeFactor / 2;
                iter_for_stepsize_factor = 0;
            }
            double numerator = stepsizeFactor * (bestObjUpBound - LBNoLift);
            double denominator = 0;
            for (LogiTask logiTask : logiTaskSet) {
                denominator = denominator
                        + ((logiTask.getTad() - logiTask.getTss()) * (logiTask.getTad() - logiTask.getTss()));
            }
            stepsize = numerator / denominator;
        } else {

            stepsize = numerator_for_stepsize / iter;
        }
////		计算stepSize的上界，确保对应各产品的所有新的lambda求和不会超过1
//		double[] sum1 = new double[para.numProd];
//		double[] sum2 = new double[para.numProd];
//		for (int p = 0; p < para.numProd; p++) {
//			sum1[p] = 0;
//			sum2[p] = para.objValFactor;
//			stepsizeUpperBound_Prod[p] = Float.MAX_VALUE;
//		}
//		for (LogiTask logiTask : logiTaskSet) {
//			sum1[logiTask.getProd().getId()] = sum1[logiTask.getProd().getId()] + logiTask.getTad() - logiTask.getTss();
//			sum2[logiTask.getProd().getId()] = sum2[logiTask.getProd().getId()] - logiTask.getLambda();
//		}
//		for (int p = 0; p < para.numProd; p++) {
//			if (sum1[p] > 0) {
//				stepsizeUpperBound_Prod[p] = sum2[p] / sum1[p];
//				stepsize = Math.min(stepsizeUpperBound_Prod[p], stepsize);
////				System.out.println("stepsizeUpperBound_Prod=(1-sum[lambda])/(sum[tad-tss])"+sum2[p]+"/"+sum1[p]+"="+stepsizeUpperBound_Prod[p]);
//			}
//		}
    }

    public void updateLambda() {
//		update lagrangian multiplers with current subgradient and stepsize  
        for (LogiTask logiTask : logiTaskSet) {
            logiTask.setLambda(Math.max(logiTask.getLambda() + stepsize * logiTask.getSubgrad(), 0));
        }
    }

    public double getObj(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq) {

        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur;
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
        double sum_taf = 0;
        double[] taf=new double[para.numProd];
        key = prodSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(para.moldChangeTime, PMKtoLogiTask.get(key).tadHeur);// 第1个产品第1个机器第1个批次的物料的开始时间
        // 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
//			timeCur=timeCur+(para.jobshopLength/para.numMach)/para.lineSpeed;
            taf[prodSeq[p]] = timeCur;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }

//        maxTaf=Math.min(Arrays.stream(taf).max().getAsDouble(),maxTaf);
//        for (int p = 0; p < para.numProd ; p++) {
//            if(tafLB[p]>taf[prodSeq[p]]){
//                System.out.println("p"+p+", tafLB="+tafLB[p]+", taf="+taf[prodSeq[p]]);
//                System.out.print("");
//            }
//        }
        return sum_taf;
    }

    public void print3() {
        for (int i = 0; i < bestLogiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < bestLogiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).tadHeur = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).tadHeur;
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
        double sum_taf = 0;
        key = bestProdSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(para.moldChangeTime, PMKtoLogiTask.get(key).tadHeur);// 第1个产品第1个机器第1个批次的物料的开始时间
        // 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < bestProdSeq.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(bestProdSeq[p]).getCycleTime() * (para.demands[bestProdSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[bestProdSeq[p]][m]; k++) {// batch从0开始编号
                    key = bestProdSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
//			timeCur=timeCur+(para.jobshopLength/para.numMach)/para.lineSpeed;
            taf[bestProdSeq[p]] = timeCur;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[bestProdSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
        System.out.println("production sequence (Heur):");
        for (int i = 0; i < bestProdSeq.length; i++) {
            System.out.print(" -> " + bestProdSeq[i] + "(" + taf[bestProdSeq[i]] + ", P"
                    + para.prodSet.get(bestProdSeq[i]).getId() + ")");
        }
        System.out.println();

        System.out.println("logistics sequence (Heur): ");
        for (int i = 0; i < bestLogiSeq.size(); i++) {
            for (int j = 0; j < bestLogiSeq.get(i).size(); j++) {
                System.out.print(" -> " + bestLogiSeq.get(i).get(j) + "(" + tad[bestLogiSeq.get(i).get(j)] + ",P"
                        + para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getProd().getId() + ",M"
                        + para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
        }

    }

    public void print(int[] bestProdSeq2, ArrayList<ArrayList<Integer>> bestLogiSeq2) {

        double[] tad = new double[para.logiTaskSet.size()];
//		double[] tss=new double[para.logiTaskSet.size()];
        double[] taf = new double[para.numProd];
        for (int i = 0; i < bestLogiSeq2.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {// 第i个路径上的第j个位置
                tad[bestLogiSeq2.get(i).get(j)] = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = tad[bestLogiSeq2.get(i).get(j)];
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
//        key = bestProdSeq2[0] + "_" + 0 + "_" + 0;
        timeCur = para.moldChangeTime;// 第1个产品第1个机器第1个批次的物料的开始时间 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < bestProdSeq2.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
//                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                    timeCur = timeCur
                            - para.prodSet.get(bestProdSeq2[p]).getCycleTime() * (para.demands[bestProdSeq2[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;
                }
                for (int k = 0; k < para.K[bestProdSeq2[p]][m]; k++) {// batch从0开始编号
                    key = bestProdSeq2[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, tad[PMKtoLogiTask.get(key).getId()]); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            taf[bestProdSeq2[p]] = timeCur;
        }

        System.out.println("production sequence (Heur):");
        for (int i = 0; i < bestProdSeq2.length; i++) {
            System.out.print(" -> " + bestProdSeq2[i] + "(" + taf[bestProdSeq2[i]] + ", P"
                    + para.prodSet.get(bestProdSeq2[i]).getId() + ")");
        }
        System.out.println();

        System.out.println("logistics sequence (Heur): ");
        for (int i = 0; i < bestLogiSeq2.size(); i++) {
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.print(" -> " + bestLogiSeq2.get(i).get(j) + "(" + tad[bestLogiSeq2.get(i).get(j)] + ",P"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getProd().getId() + ",M"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
        }

    }

    public void print2(int[] bestProdSeq2, ArrayList<ArrayList<Integer>> bestLogiSeq2) {// 输出tss和tad，看哪里的充电可以提前
        double[] tad = new double[para.logiTaskSet.size()];
        double[] tss = new double[para.logiTaskSet.size()];
        double[] tssNoLogi = new double[para.logiTaskSet.size()];
        double[] tssLate = new double[para.logiTaskSet.size()];
        double[] taf = new double[para.numProd];
        for (int i = 0; i < bestLogiSeq2.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {// 第i个路径上的第j个位置
                tad[bestLogiSeq2.get(i).get(j)] = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = tad[bestLogiSeq2.get(i).get(j)];
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
        key = bestProdSeq2[0] + "_" + 0 + "_" + 0;
        timeCur = tad[PMKtoLogiTask.get(key).getId()];// 第1个产品第1个机器第1个批次的物料的开始时间 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < bestProdSeq2.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur
                            - para.prodSet.get(bestProdSeq2[p]).getCycleTime() * (para.demands[bestProdSeq2[p]] - 1);
                }
                for (int k = 0; k < para.K[bestProdSeq2[p]][m]; k++) {// batch从0开始编号
                    key = bestProdSeq2[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, tad[PMKtoLogiTask.get(key).getId()]); // 真正的开始时间是物料到达和最早开工时间取大；
                    tss[PMKtoLogiTask.get(key).getId()] = timeCur;
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            taf[bestProdSeq2[p]] = timeCur;
        }

        double[][] tss_PM = new double[para.numProd][para.numMach];// 产品p在机器m上的基准时间
        double[] tss_P = new double[para.numProd];// 产品p的基准时间
        tss_P[bestProdSeq2[0]] = para.getDistance(para.depot, para.machSet.get(0)) * 2 * para.consumingRate
                * para.chargingRate + para.getDistance(para.depot, para.machSet.get(0)) / para.vehSpeed;
        para.prodSet.get(bestProdSeq2[0]).setTaf(tss_P[bestProdSeq2[0]]
                + (para.demands[bestProdSeq2[0]] + para.numMach - 1) * para.cycTime[bestProdSeq2[0]]);
        for (int i = 1; i < para.numProd; i++) { // 下一个产品的开始时间由上一个产品的开始时间决定
            tss_P[bestProdSeq2[i]] = tss_P[bestProdSeq2[i - 1]]
                    + (para.demands[bestProdSeq2[i - 1]] + para.numMach - 1) * para.cycTime[bestProdSeq2[i - 1]];
            para.prodSet.get(bestProdSeq2[i]).setTaf(tss_P[bestProdSeq2[i]]
                    + (para.demands[bestProdSeq2[i]] + para.numMach - 1) * para.cycTime[bestProdSeq2[i]]);
        }
        for (int i = 0; i < para.numProd; i++) {
            for (int j = 0; j < para.numMach; j++) {// machine从0开始编号
                tss_PM[bestProdSeq2[i]][j] = tss_P[bestProdSeq2[i]] + j * para.cycTime[bestProdSeq2[i]];
            }
        }

        // 更新tss
        double tss_PMK = 0;
        for (int i = 0; i < para.numProd; i++) {// product从0开始编号
            for (int j = 0; j < para.numMach; j++) {// machine从0开始编号
                for (int k = 0; k < para.K[bestProdSeq2[i]][j]; k++) {// batch从0开始编号
                    tss_PMK = tss_PM[bestProdSeq2[i]][j]
                            + k * para.quanProdMach[bestProdSeq2[i]][j] * para.cycTime[bestProdSeq2[i]];
                    tssNoLogi[PMKtoLogiTask.get(bestProdSeq2[i] + "_" + j + "_" + k).getId()] = tss_PMK;
                    tssLate[PMKtoLogiTask.get(bestProdSeq2[i] + "_" + j + "_" + k).getId()] = taf[bestProdSeq2[i]]
                            + ((k + 1 - 1) * para.quanProdMach[bestProdSeq2[i]][j] - para.numMach + (j + 1)
                            - para.demands[bestProdSeq2[i]]) * para.cycTime[bestProdSeq2[i]];
                }
            }
        }

        System.out.println("production sequence:");
        for (int i = 0; i < bestProdSeq2.length; i++) {
            System.out.print(" -> " + bestProdSeq2[i] + "(" + taf[bestProdSeq2[i]] + ", P"
                    + para.prodSet.get(bestProdSeq2[i]).getId() + ")");
        }
        System.out.println();

        System.out.println("logistics sequence: ");
        for (int i = 0; i < bestLogiSeq2.size(); i++) {
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.print(" -> " + bestLogiSeq2.get(i).get(j) + "(" + tad[bestLogiSeq2.get(i).get(j)] + ","
                        + tss[bestLogiSeq2.get(i).get(j)] + ",P"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getProd().getId() + "M"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
        }

        // 生产信息输出
        System.out.printf("%-16s", "ProdID");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10d", bestProdSeq2[p]);
        }
        System.out.println();
        System.out.printf("%-16s", "taf");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10.2f", taf[bestProdSeq2[p]]);
        }
        System.out.println();
        System.out.println();

        // 物流信息输出

        for (int i = 0; i < bestLogiSeq2.size(); i++) {
            System.out.printf("%-16s", "LogiID");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10s",
                        bestLogiSeq2.get(i).get(j) + "(P"
                                + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getProd().getId() + "M"
                                + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
            System.out.printf("%-16s", "tad");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tad[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.printf("%-16s", "tss");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tss[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.printf("%-16s", "tssNoLogi");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tssNoLogi[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.printf("%-16s", "tssLate");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tssLate[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.println();
        }

    }
}
