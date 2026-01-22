package column_generation;

import ilog.concert.IloException;
import ilog.cplex.IloCplex;

public class Parameters {

public static void configureCplex(RstMasterProblem masterproblem) {
	try {
		// branch and bound
		masterproblem.rmpModel.setParam(IloCplex.Param.MIP.Strategy.NodeSelect, 1);
		masterproblem.rmpModel.setParam(IloCplex.Param.MIP.Strategy.Branch,1);
		//masterproblem.cplex.setParam(IloCplex.Param.Preprocessing.Presolve, true);
		// display options
		masterproblem.rmpModel.setParam(IloCplex.Param.MIP.Display, 2);
		masterproblem.rmpModel.setParam(IloCplex.Param.Tune.Display, 1);
		masterproblem.rmpModel.setParam(IloCplex.Param.Simplex.Display, 0);
	}
	catch (IloException e) {System.err.println("Concert exception caught: " + e);}
}

}