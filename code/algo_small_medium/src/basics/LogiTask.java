package basics;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

public class LogiTask implements Serializable{


	private static final long serialVersionUID = 2839132336012975496L;
	private int id;


	private Machine correspMach;
	private Product correspProd;
	private int batchId;
	private double transTime;
	private double lambda; //对应的是拉格朗日乘子
	private double dualVal;//对应列生成的对偶变量
	private double tss;//tssEarl_no_dist
	public double tss_with_dist_Earl;
	public double tss_with_dist_Late;
	public double tssHeur_Late;
	public double tssHeur_Earl;
	public double tssHeur;
	public double tadHeur;
	private double tad;
	private double subgrad;
	private int batchQuant;
	private double assemTime;
	private double processTime;//对应列生成算法中的travelTime+chargingTime
	private double chargingTime;
	
    public double getChargingTime() {
		return chargingTime;
	}

	public void setChargingTime(double chargingTime) {
		this.chargingTime = chargingTime;
	}

	@Override
    public String toString() {
        return "[id" + this.id  + ": for Product " + this.correspProd.getId() 
        + ", Machine "+this.correspMach.getId()+", Batch "+this.batchId
        + ", transport time " + this.getTransTime()+", lagMultiplier "+this.getLambda()+ ", ratio " + this.lambda/this.getTransTime()+"]";
    }
    
    //序列化与反序列化的方式实现深度拷贝
    public Object deepCopy() throws IOException, ClassNotFoundException{
        //字节数组输出流，暂存到内存中
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        //序列化
        ObjectOutputStream oos = new ObjectOutputStream(bos);
        oos.writeObject(this);
        ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        ObjectInputStream ois = new ObjectInputStream(bis);
        //反序列化
        return ois.readObject();
    }

	public void setId(int id) {
		this.id=id;
	}
	
	public void setBatchId(int batchId) {
		this.batchId=batchId;
	}
	
	public int getId() {
		return id;
	}

	public int getBatchId() {
		return batchId;
	}
	
	public void setProd(Product prod) {
		this.correspProd=prod;
	}
	
	public Product getProd() {
		return this.correspProd;
	}
	
	public void setMach(Machine mach) {
		this.correspMach=mach;
	}
	
	public Machine getMach() {
		return this.correspMach;
	}

	public void setDualValue(double dualVal) {
		this.dualVal=dualVal;
	}
	
	public double getDualValue() {
		return this.dualVal;
	}

	public double getTss() {
		return tss;
	}

	public void setTss(double tss) {
		this.tss = tss;
	}

	public double getTad() {
		return tad;
	}

	public void setTad(double tad) {
		this.tad = tad;
	}

	public int getBatchQuant() {
		return batchQuant;
	}

	public void setBatchQuant(int batchQuant) {
		this.batchQuant = batchQuant;
	}

	public double getSubgrad() {
		return subgrad;
	}

	public void setSubgrad(double subgrad) {
		this.subgrad = subgrad;
	}

	public double getAssemTime() {
		return assemTime;
	}

	public void setAssemTime(double assemTime) {
		this.assemTime = assemTime;
	}

	public double getTransTime() {
		return transTime;
	}

	public void setTransTime(double transTime) {
		this.transTime = transTime;
	}

	public double getLambda() {
		return lambda;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	public void setProcessingTime(double processTime) {
		this.processTime = processTime;
	}
	
	public double getProcessTime() {
		return processTime;
	}


}
