package test;

import java.util.ArrayList;
import java.util.Random;

public class TestTableOutput {
	public double ini_heur(int[] prodSeq, ArrayList<ArrayList<ArrayList<Integer>>> logiSeq_rand) {
		double bestObjVal = Float.MAX_VALUE;
		for (int ii = 0; ii < 2000; ii++) {
			if (Math.random() < 0.5) {
				Random rnd = new Random();
				int p_change = rnd.nextInt(para.numProd - 1);
				int temp = prodSeq[p_change + 1];
				prodSeq[p_change + 1] = prodSeq[p_change];
				prodSeq[p_change] = temp;
			}
			ArrayList<ArrayList<ArrayList<Integer>>> cumb_logiSeq_rand = new ArrayList<ArrayList<ArrayList<Integer>>>();
			for (int n = 0; n < num_sample; n++) {
				ArrayList<ArrayList<Integer>> cumb_logiSeq = new ArrayList<ArrayList<Integer>>();
				for (int r = 0; r < para.numVeh; r++) {
					cumb_logiSeq.add(new ArrayList<Integer>());
				}
				cumb_logiSeq_rand.add(cumb_logiSeq);
			}

			ArrayList<ArrayList<Integer>> first_rand = new ArrayList<ArrayList<Integer>>();// prodSeq中对每个订单p对应的第一个物流任务编号
			ArrayList<ArrayList<Integer>> last_rand = new ArrayList<ArrayList<Integer>>();
			for (int n = 0; n < num_sample; n++) {
				ArrayList<Integer> first = new ArrayList<>();// prodSeq中对每个订单p对应的第一个物流任务编号
				ArrayList<Integer> last = new ArrayList<>();
				for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号3.
					int pindex = prodSeq[p];
					ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI_sample.get(n).get(pindex);
					ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
					for (int index = 0; index < arrayofPMKI.size(); index++) {
						arrayOfI.add(arrayofPMKI.get(index)[3]);
					}
					first.add(arrayOfI.get(0));
					last.add(arrayOfI.get(arrayOfI.size() - 1));
				}
				first_rand.add(first);
				last_rand.add(last);
			}

			ArrayList<Integer> prodPriority = new ArrayList<Integer>();
			for (int p = 0; p < prodSeq.length; p++) {
				prodPriority.add(prodSeq[p]);
			}
			int[][] priorityLogiTask = new int[maxLogiTaskSize][num_sample];// 每个场景中logiTask的优先级,按生产顺序排的
			for (int n = 0; n < num_sample; n++) {
				for (int i = 0; i < logiSeq_rand.get(n).size(); i++) {
					for (int j = 0; j < logiSeq_rand.get(n).get(i).size(); j++) {
						int logiTaskId = logiSeq_rand.get(n).get(i).get(j);
						int prodId = para.ItoLogiTask_sample.get(n).get(logiTaskId).getProd().getId();
						priorityLogiTask[logiTaskId][n] = prodPriority.indexOf(prodId);
					}
				}
			}

			for (int n = 0; n < num_sample; n++) {
				for (int i = 0; i < logiSeq_rand.get(n).size(); i++) {
					ArrayList<LogiTaskPrior> priorSet = new ArrayList<LogiTaskPrior>();
					for (int j = 0; j < logiSeq_rand.get(n).get(i).size(); j++) {
						int logiTaskId = logiSeq_rand.get(n).get(i).get(j);
						priorSet.add(new LogiTaskPrior(logiTaskId, priorityLogiTask[logiTaskId][n]));
					}
					Collections.sort(priorSet, new LogiTaskSortByPrior());
					for (int j = 0; j < priorSet.size(); j++) {
						cumb_logiSeq_rand.get(n).get(i).add(priorSet.get(j).id);
					}
				}
			}

			bestObjVal = Math.min(getObj_rand(prodSeq, cumb_logiSeq_rand), getObj_rand(prodSeq, logiSeq_rand));
			ArrayList<ArrayList<ArrayList<Integer>>> logiSeq_rand_new = listCopy3(cumb_logiSeq_rand);

			for (int n = 0; n < num_sample; n++) {
				int iter = 0;
				Random rnd = new Random();
				while (iter < 20) {
					int ov;
					int op;
					while (true) {
						int a = rnd.nextInt(para.numVeh);
						if (logiSeq_rand_new.get(n).get(a).size() > 1) {
							ov = a;
							op = rnd.nextInt(logiSeq_rand_new.get(n).get(a).size());
							break;
						}
					}
					int o_LogiTaskId = logiSeq_rand_new.get(n).get(ov).get(op);
					int o_ProdId = para.ItoLogiTask_sample.get(n).get(o_LogiTaskId).getProd().getId();
					int o_Priority = prodPriority.indexOf(o_ProdId);

					HashSet<int[]> feasInnerInsertPos = new HashSet<>();// 相同优先级的插入（插入在指定位置之前）
					HashSet<int[]> feasAfterInsertPos = new HashSet<>();// 插入在指定位置之后
					HashSet<int[]> feasEndInsertPos = new HashSet<>();// 插入在队列末尾
					HashSet<int[]> feasSwapPos = new HashSet<>();// 相同优先级的交换

					for (int i = 0; i < logiSeq_rand_new.get(n).size(); i++) {
						for (int j = 0; j < logiSeq_rand_new.get(n).get(i).size(); j++) {
							int[] fea = { i, j };
							int n_LogiTaskId = logiSeq_rand_new.get(n).get(i).get(j);
							if (o_LogiTaskId != n_LogiTaskId) {
								int n_ProdId = para.ItoLogiTask_sample.get(n).get(n_LogiTaskId).getProd().getId();
								int n_Priority = prodPriority.indexOf(n_ProdId);
								if (n_Priority == o_Priority) {
									feasSwapPos.add(fea);
									feasInnerInsertPos.add(fea);
								}
								if (j < logiSeq_rand_new.get(n).get(i).size() - 1) {
									int n2_LogiTaskId = logiSeq_rand_new.get(n).get(i).get(j + 1);
									int n2_ProdId = para.ItoLogiTask_sample.get(n).get(n2_LogiTaskId).getProd().getId();
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
					while (!improved && (!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty()
							|| !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty())) {
						switch (rndOperator.nextInt(4)) {
						case 0:
							if (!feasInnerInsertPos.isEmpty()) {
								Iterator<int[]> iterator = feasInnerInsertPos.iterator();
								int[] fea = iterator.next();
								feasInnerInsertPos.remove(fea);
								ArrayList<ArrayList<ArrayList<Integer>>> neigh_logiSeq = listCopy3(logiSeq_rand_new);
								innerInsertOperator(neigh_logiSeq.get(n), ov, op, fea[0], fea[1]);
								double testObjVal = getObj_rand(prodSeq, neigh_logiSeq);
								if (testObjVal <= bestObjVal) {
									bestObjVal = testObjVal;
									logiSeq_rand_new = neigh_logiSeq;
									improved = true;
								}
							}
							break;
						case 1:
							if (!feasEndInsertPos.isEmpty()) {
								Iterator<int[]> iterator = feasEndInsertPos.iterator();
								int[] fea = iterator.next();
								feasEndInsertPos.remove(fea);
								ArrayList<ArrayList<ArrayList<Integer>>> neigh_logiSeq = listCopy3(logiSeq_rand_new);
								endInsertOperator(neigh_logiSeq.get(n), ov, op, fea[0], fea[1]);
								double testObjVal = getObj_rand(prodSeq, neigh_logiSeq);
								if (testObjVal <= bestObjVal) {
									bestObjVal = testObjVal;
									logiSeq_rand_new = neigh_logiSeq;
									improved = true;
								}
							}
							break;
						case 2:
							if (!feasAfterInsertPos.isEmpty()) {
								Iterator<int[]> iterator = feasAfterInsertPos.iterator();
								int[] fea = iterator.next();
								feasAfterInsertPos.remove(fea);
								ArrayList<ArrayList<ArrayList<Integer>>> neigh_logiSeq = listCopy3(logiSeq_rand_new);
								afterInsertOperator(neigh_logiSeq.get(n), ov, op, fea[0], fea[1]);
								double testObjVal = getObj_rand(prodSeq, neigh_logiSeq);
								if (testObjVal <= bestObjVal) {
									bestObjVal = testObjVal;
									logiSeq_rand_new = neigh_logiSeq;
									improved = true;
								}
							}
							break;
						case 3:
							if (!feasSwapPos.isEmpty()) {
								Iterator<int[]> iterator = feasSwapPos.iterator();
								int[] fea = iterator.next();
								feasSwapPos.remove(fea);
								ArrayList<ArrayList<ArrayList<Integer>>> neigh_logiSeq = listCopy3(logiSeq_rand_new);
								swapOperator(neigh_logiSeq.get(n), ov, op, fea[0], fea[1]);
								double testObjVal = getObj_rand(prodSeq, neigh_logiSeq);
								if (testObjVal <= bestObjVal) {
									bestObjVal = testObjVal;
									logiSeq_rand_new = neigh_logiSeq;
									improved = true;
								}
							}
							break;
						}
					}
					iter++;
				}

				int[][] posi = new int[para.numProd][4];
				for (int i = 0; i < logiSeq_rand_new.get(n).size(); i++) {
					for (int j = 0; j < logiSeq_rand_new.get(n).get(i).size(); j++) {
						int logiTaskId = logiSeq_rand_new.get(n).get(i).get(j);
						int prodId = para.ItoLogiTask_sample.get(n).get(logiTaskId).getProd().getId();
						if (first_rand.get(n).contains(logiTaskId)) {
							posi[prodId][0] = i;
							posi[prodId][1] = j;
						}
						if (last_rand.get(n).contains(logiTaskId)) {
							posi[prodId][2] = i;
							posi[prodId][3] = j;
						}
					}
				}
				ArrayList<ArrayList<ArrayList<Integer>>> neigh_logiSeq = listCopy3(logiSeq_rand_new);
				bestObjVal = Math.min(bestObjVal, LocalImprove(prodSeq, neigh_logiSeq, bestObjVal, posi, n));
			}
		}

		return bestObjVal;
	}

//	public static void main(String[] args) {
//		// TODO Auto-generated method stub
////		System.out.println("星期一\t星期二\t星期三\t星期四\t星期五\t星期六\t星期日");
////		for (int j = 1; j <= 30; j++) {
////			// 用两个制表符调整位置
////			System.out.printf(1234567 + "\t");
////			if (j % 7 == 0) {
////				System.out.println();
////			}
////		}
//		System.out.printf("%-11s\t","TaskID");
//		for(int i=0;i<6;i++) {
//			System.out.printf("%10d",i);
//		}
//		System.out.println();
//		
//		System.out.printf("%-11s\t","lambda");
//		for(int i=0;i<6;i++) {
//			System.out.printf("%10.2f",IdToLogiTask);
//		}
//		System.out.println();
//	}
	
	public int a(int s) {
		if(s==1) {
			System.out.println("-100");
			return -100;
		}
		else {
			System.out.println("100");
			return 100;
		}
	}

}
