import java.util.ArrayList;

/** Moving window interval to calculate cMBFs from
 */
public class IntWindow {
	//-- Static parameters --
	private static int intSize;
	private static double medianMult;
	
	/** Window size in terms of lines **/
	private int windowSize;
	//-- Window stats/parameters -- 
	private String chromNum;
	private int absStart;
	private int start;
	private int pos;
	private int end;
	private int relIndex;
	private boolean endOfChrom;

	/** ArrayList of IntStats in the window */
	protected ArrayList<IntStats> winRA;

	/** Constructor
	 * @param wS - window size in terms of lines
	 * @param iS - interval size
	 */
	public IntWindow(String cN, int wS, int iS, double mM, int p) {
		chromNum = cN;
		intSize = iS;
		windowSize = wS / intSize;
		winRA = new ArrayList<IntStats>(windowSize);
		
		medianMult = mM;
		endOfChrom = false;
		pos = p;
		
		absStart = p;
		start = p;
		relIndex = 0;
	}
	
	//-- Get window statistics -- 
	public String getChromNum() {
		return chromNum;
	}
	
	public int getIndexStart() {
		return winRA.get(relIndex).getStart();
	}

	public int getIndexEnd() {
		return winRA.get(relIndex).getEnd();
	}

	public double getIndexReadCount() {
		return winRA.get(relIndex).getReadCount();
	}

	public double getIndexcMBF() {
		return winRA.get(relIndex).getcMBF();
	}
	
	public int getLastStartIndex() {
		return winRA.get(winRA.size() - 1).getStart();
	}
	
	/** Set end of chromosome */
	public void setEndOfChrom() {
		endOfChrom = true;
	}

	/** Calculate expected start and end */
	private void calcStartEnd() {
		if (!endOfChrom) { //don't change window start/end if reached end of file
			if (windowSize % 2 != 0) { //if window size odd
				start = pos - intSize * (windowSize/ 2);
				end = pos + intSize * ((windowSize + 1) / 2);
			}
			else {
				start = (pos - (int) (Math.floor((windowSize - 1)/2) * intSize));
				end = (pos + (int) (Math.ceil((windowSize + 1)/2) * intSize));
			}
		}
		
		if (start < absStart) {
			end = (end - start) + absStart;
			start = absStart;
		}
		
		relIndex = (pos - start) / intSize;
	}
	
	/** Used for when window is not filled but reached end of file - assumingly 0 - [last line's start] */
	public int setSmallerWindowSize() {
		windowSize = winRA.size();
		end = winRA.get(winRA.size() - 1).getEnd();
		return windowSize;
	}

	/**
	 * Calculates the number of elements left to fill in window
	 * @return
	 */
	public int toFill() {
		return windowSize - winRA.size();		
	}

	public boolean full() {
		return winRA.size() == windowSize;
	}

	/** Insert an IntStats object to window */
	public boolean insert(IntStats is) {
		if (winRA.size() >= windowSize) { //full
			throw new IndexOutOfBoundsException();
		}
		return winRA.add(is);
	}

	//Available for changes in calculation
	@SuppressWarnings("unused")
	private double calcMean() { //used in place of Z
		double sum = 0;
		int count = 0;
		for (int i = 0; i < windowSize; i++) {
			sum += winRA.get(i).getReadCount();
			count++;
		}
		return sum / count;
	}

	private double calcMedian() {
		@SuppressWarnings("unchecked")
		ArrayList<IntStats> sortedRA = (ArrayList<IntStats>) winRA.clone();
		sortedRA.sort(new ISWGSepReadCountSorter());
		int middle = 0;
		if (windowSize % 2 == 0) { //if even, average middle two entries
			middle = windowSize/2 - 1;
			return (sortedRA.get(middle).getReadCount() + sortedRA.get(middle + 1).getReadCount()) / 2;
		}
		else {
			middle = windowSize/2;
			return sortedRA.get(middle).getReadCount(); //else odd, return middle entry
		}
		
	}

	/** Calculate the cMBF (1 - exp(-Z^2/2), Z = RC/M ) for the position using the median of the window */
	public double calccMBF() throws Exception{ //currently uses mean
		double noiseEst = medianMult * calcMedian();
		if (noiseEst == 0) {
			throw new Exception("Noise level estimated to be 0: \n"
					+ chromNum + "\t" + start + "\t" + end + "\t" + pos);
		}
		
		double Z = getIndexReadCount() / noiseEst;
		double cmbf = 1 - Math.exp(-1 * (Z*Z) / 2);
		winRA.get(relIndex).setcMBF(cmbf);
		return cmbf;
	}

	/** Increment center and move it along the chromosome */
	public void incrCenter() {
		pos += intSize;
		calcStartEnd();
		if (!endOfChrom && winRA.get(0).getStart() < start) { //remove first element in RA if moving out of range and haven't reached end of file
			winRA.remove(0);
		}
	}

	//Used for debugging
	@Override
	public String toString() {
		return "windowSize:" + windowSize + ", pos:" + pos + ", start:" + start + ", end:" + end + ", relIndex" + relIndex + "\nwinRA:" + winRA.toString();
	}

}
