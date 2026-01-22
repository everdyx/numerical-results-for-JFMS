package main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import basics.Instance;
import subgradient_method.ProdRatio;
import subgradient_method.TssWithoutLogiCons;

/**
 * Sol
 */
public class IntegratedSol {
    public int[] prodSeq;
    public ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
    // public ArrayList<ArrayList<Double>> charging=new ArrayList<ArrayList<Double>>();
    public double[] charging;
    public int[] quantity;
    
    int numProd;
    int numVeh;
    Instance para;

    public IntegratedSol(int numProd, int numVeh, Instance para){
        this.numProd=numProd;
        this.numVeh=numVeh;
        this.prodSeq=new int[numProd];
        for (int i = 0; i < numVeh; i++) {
            this.logiSeq.add(new ArrayList<Integer>());
            // this.charging.add(new ArrayList<Double>());
        }
        this.charging=new double[para.logiTaskSet.size()];
        this.quantity=new int[para.logiTaskSet.size()];
        this.para=para;
    }

    /**
     * 
     */
    public void inilizatialize(int[] prodSeq_fcfs,ArrayList<ArrayList<Integer>> logiSeq_fcfs){
        Random rnd=new Random();

        // for (int i = 0; i < prodSeq.length; i++) {
        //     prodSeq[i]=i;
        // }
        
        // for (int i = 0; i < para.logiTaskSet.size(); i++) {
        //     int index=rnd.nextInt(numVeh);
        //     // System.out.println(para.logiTaskSet.get(i).getId());
        //     logiSeq.get(index).add(para.logiTaskSet.get(i).getId());
        // }
        prodSeq=prodSeq_fcfs;
        for (int i = 0; i < logiSeq_fcfs.size(); i++) {
            for (int j = 0; j < logiSeq_fcfs.get(i).size(); j++) {
                logiSeq.get(i).add(logiSeq_fcfs.get(i).get(j));
            }
        }
        for (int i = 0; i < charging.length; i++) {
            charging[i]=rnd.nextDouble(0,20*60+1);
        }
        for (int p = 0; p < para.numProd; p++) {
            for (int m = 0; m < para.numMach; m++) {
                int accu=0;
                int ub=para.quanProdMach[p][m];
                for (int k = 0; k< para.K[p][m]; k++) {
                    int lid=para.PMKtoLogiTask.get(p+"_"+m+"_"+k).getId();
                    if (para.K[p][m]==1) {
                        quantity[lid]=para.demands[p];    
                    }else{
                        if (k==para.K[p][m]-1) {
                            quantity[lid]=para.demands[p]-accu;
                        }else {
                            quantity[lid]=Math.min(rnd.nextInt(para.demands[p]+1),ub);
                            accu+=quantity[lid];
                        }
                    }
                }
            }
        }
    }
    public IntegratedSol copy(){
        IntegratedSol sol2=new IntegratedSol(this.prodSeq.length, this.logiSeq.size(), this.para);
        for (int i = 0; i < this.prodSeq.length; i++) {
            sol2.prodSeq[i] = this.prodSeq[i];
        }
        
        for (int i = 0; i < this.logiSeq.size(); i++) {
            for (int j = 0; j < this.logiSeq.get(i).size(); j++) {
                sol2.logiSeq.get(i).add(this.logiSeq.get(i).get(j));
            }
        }
        for (int i = 0; i < this.charging.length; i++) {
            sol2.charging[i] = this.charging[i];
        }
        // for (int i = 0; i < this.logiSeq.size(); i++) {
        //     for (int j = 0; j < this.logiSeq.get(i).size(); j++) {
        //         sol2.charging.get(i).add(this.charging.get(i).get(j));
        //     }
        // }
        for (int i = 0; i < this.quantity.length; i++) {
            sol2.quantity[i] = this.quantity[i];
        }
        return sol2;
    }

    public void printSol(){
        System.out.println("Prod Sequence:");
        for (int i = 0; i < prodSeq.length; i++) {
            System.out.print(prodSeq[i]+" ");
        }
        System.out.println();
        System.out.println("Logi Sequence:");
        for (int i = 0; i < logiSeq.size(); i++) {
            for (int j = 0; j < logiSeq.get(i).size(); j++) {
                System.out.print(logiSeq.get(i).get(j)+" ");
            }
            System.out.println();
        }
    }
}