package dynamic_programming;

public class StageState {
	private int stage;
	private int state;

	public StageState() {
		// TODO Auto-generated constructor stub
	}
	
	public StageState(int stage, int state) {
		this.stage=stage;
		this.state=state;
	}
	
	public int getStage() {
		return stage;
	}

	public void setStage(int stage) {
		this.stage = stage;
	}

	public int getState() {
		return state;
	}

	public void setState(int state) {
		this.state = state;
	}

	public String toString() {
		return stage + "-" + state;
	}
    
	@Override
    public int hashCode()
    {
        return this.toString().hashCode();
    }
	

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || this.getClass() != o.getClass()) {
			return false;
		}

		StageState stage_state = (StageState) o;

		if (this.stage==stage_state.stage && this.state==stage_state.state) { 
			return true;
		}
		else {
			return true;
		}
	}


}
