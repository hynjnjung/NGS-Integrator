 package nih.nhlbi.esbl.ngs;


import java.util.Comparator;

/**
 * Object that holds the read count and cMBF of a position in the chromosome
 * Can get/set readCount/cMBF
 */
class IntStats {
	private static double defZero = 0.1;
	
	private String _chromNum;
	private int _start;
	private int _end;
	private int _width;
	private double _readCount;
	private double _cMBF;
	
	/** Constructor */
	public IntStats(String cN, int s, int e, int w, int rc) {
		_chromNum = cN;
		_start = s;
		_end = e;
		_width = w;

		setReadCount(rc);
		
		_cMBF = -1;
	}
	
	//-- Set/get parameters --
	public static void setdefZero(double dz) {
		defZero = dz;
	}
	
	public String getChromNum() {
		return _chromNum;
	}
	
	public int getStart() {
		return _start;
	}
	
	public int getWidth() {
		return _width;
	}
	
	public int getEnd() {
		return _end;
	}
	
	public double getReadCount() {
		return _readCount;
	}
	
	public double getcMBF() {
		return _cMBF;
	}
	
	public void setReadCount(int rc) {
		if (rc > 0) {
			_readCount = rc;
		}
		else if (rc == 0){
			_readCount = defZero;
		}
		else {
			System.err.println("readCount cannot be negative" + rc);
			System.exit(1);
		}
	}
	
	
	public void setcMBF(double cmbf) throws Exception {
		if (cmbf >= 0 && cmbf <= 1) {
			_cMBF = cmbf;
		}
		else{
			throw new Exception("cMBF cannot be outside the range of 0 to 1: " + cmbf + "\n"
					+ this._chromNum + "\t" + this._start + "\t" + this._end + "\t" + this._readCount);
		}
	}
	
	@Override
	public String toString() {
		return "[start:" + _start + ", end:" + _end + ", readCount:" + _readCount + ", cMBF:" + _cMBF + "]";
	}
	
}

/** Comparator class for sorting the IntStats in order to find the median*/
class ISWGSepReadCountSorter implements Comparator<IntStats> {

	@Override
	public int compare(IntStats ps1, IntStats ps2) {
		return (int) (10*(ps1.getReadCount() - ps2.getReadCount()));
	}
	
}