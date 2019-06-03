package nih.nhlbi.esbl.ngs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/** Calculates the cMBF for a given input genome coverage file.*/
public class Calculator {
	/** Window size for calculating cMBF (in bp) */
	private static Integer windowbpSize;
	/** Median multiple */
	private static double medianMult;
	/** Minimum read count (in place of zero) */
	private static double defZero;
	
	private static String inBed;
	
	private static String baseName;
	private static File chromDir;
	private static String chromDirName;
	private static String outBaseName;
	private static String spacer;
	
	/**
	 * Main method, entry point for jar
	 *@param args - interval size, window size, medianMult, minRC, {input .bed file, output base} names
	 * !! Window size must be a multiple of interval size !!
	 * !! Ignores any lines that do not parse as a data line (headers, etc) !!
	 */
	public static void main(String[] args) throws Exception{
		//-- Parse options --
		parseOptions(args);
		
		//-- Split chromosome files --
		ArrayList<File> chromFiles = null;
		try {
			chromFiles = splitChromFiles(inBed);
		} catch (FileNotFoundException e) {
			System.err.println(e.getMessage());	
			System.exit(1);
		}
		
		
		//-- Parallel process chromosomes -- 
		int nThreads = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		
		ArrayList<Future<File>> futureList = new ArrayList<Future<File>>();
		ArrayList<File> outFileList = new ArrayList<File>();
		
		// Execute ProcessChromFiles for each file
		for (File f : chromFiles) {
			Callable<File> worker = new ProcessChromFile(chromDirName, f);
			Future<File> submit = executor.submit(worker);
			futureList.add(submit);
		}
		
		// Get the resulting output files from each process for concatenation
		for (Future<File> future : futureList) {
            try {
            	File f = future.get();
                System.out.println(f.getPath() + " completed");
            	outFileList.add(f);
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
		
		executor.shutdown();
		executor.awaitTermination(12, TimeUnit.HOURS);
		
		//-- Concatenate files --
		try {
			mergeWholeChrFile(outFileList);
		} catch (FileNotFoundException e) {
			System.err.println(e.getMessage());
		}
		
		//-- Clean up Chrom files and directory--
		boolean clean = true;
		for (File f: chromFiles) {
			clean = clean && f.delete();
		}
		clean = clean && chromDir.delete();
		
		if (!clean) {
			System.out.println("Unable to delete chromosome files.");
		}
		else {
			System.out.println("Split chromosome bed files cleaned.");
		}
		
		System.out.println("Completed.");
	}
	
	/**
	 * Parse commandline options using CLI
	 * @param args - commandline arguments from main method
	 */
	public static void parseOptions(String[] args) {
		//-- Set up CLI option parsing -- 
		Options options = new Options();
		
		Option infile = new Option("i", "input", true, "[req] input file path, must be a bed file");
		infile.setRequired(true);
		options.addOption(infile);
		
		Option outfileName = new Option("o", "output", true, "[opt] output file base name");
		outfileName.setRequired(false);
		options.addOption(outfileName);
		
		Option winSize = new Option("w", "windowSize", true, "[req] window size for calculating cMBF, must be a multiple of interval size");
		winSize.setRequired(true);
		options.addOption(winSize);
		
		Option medMult = new Option("m", "medMult", true, "[opt] median multiple, cannot be 0");
		medMult.setRequired(false);
		options.addOption(medMult);
		
		Option dZ = new Option("z", "defaultZero", true, "[opt] default number to replace zero");
		dZ.setRequired(false);
		options.addOption(dZ);
		
		Option help = new Option("h", "help", false, "");
		options.addOption(help);
		
		//-- Parse options --
		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd = null;
		
		String header = "Calculate cMBFs from NGS genome read count data \n\n"; 
		String footer = "\nPlease see the GitHub at https://github.com/ESBL/NGS-IT for more information.\n"
				+ "Questions or issues can be directed to ESBL at https://esbl.nhlbi.nih.gov/contact.html";
		
		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("Calculator", header, options, footer, false);
			
			System.exit(1);
		}

		if (cmd.hasOption("h")) {
			formatter.printHelp("Calculator", header, options, footer, false);
		}
		
		//Default options
		outBaseName = "out";
		windowbpSize = 10000;
		medianMult = 1;
		defZero = 0.5;
		
		//File options
		inBed = cmd.getOptionValue("i");
		if (cmd.hasOption("o")) {
			outBaseName = cmd.getOptionValue("o");
		}
		if (!inBed.endsWith(".bed")) {
			System.err.println("Incorrect input filetype (must be a .bed).");
			System.exit(0);
		}
		
		//Other parameters
		if (cmd.hasOption("w")) {
			windowbpSize = Integer.parseInt(cmd.getOptionValue("w"));
		}
		if (cmd.hasOption("m")) {
			double mm = Double.parseDouble(cmd.getOptionValue("m"));
			if (mm <= 0) {
				System.err.println("Median multiple must be greater than 0.");
				System.exit(0);
			}
			medianMult = mm;
		}
		if (cmd.hasOption("z")) {
			double mrc = Double.parseDouble(cmd.getOptionValue("z"));
			if (mrc <= 0 || mrc >= 1) {
				System.err.println("Default zero must be 0 < dZ < 1.");
				System.exit(0);
			}
			defZero = mrc;
		}
		
		//-- Set static parameters --
		IntStats.setdefZero(defZero);
		ProcessChromFile.setWinSize(windowbpSize);
		ProcessChromFile.setMedMult(medianMult);
		ProcessChromFile.setOutBaseName(outBaseName);
	}
	
	/**
	 * Given a genome-wide input .bed file, split into smaller chromosome .bed files
	 * @param infile - input .bed file
	 * @return ArrayList of the split chromosome .bed files
	 * @throws FileNotFoundException - if any files cannot be found, read or written
	 */
	public static ArrayList<File> splitChromFiles(String infile) throws FileNotFoundException {
		//-- Set up scanner for reading input file--
		ArrayList<File> chromFiles = new ArrayList<File>();
		Scanner sc;
		try {
			sc = new Scanner(new File(infile));
		} catch (FileNotFoundException e) {
			throw new FileNotFoundException("File not found: " + infile);
		}
		
		//-- Set up Chrom directory and file names--
		baseName = infile.substring(0, infile.length() - 4);
		chromDirName = baseName + "_ChromFiles";
		spacer = "/";
		
		chromDir = new File(chromDirName);
		boolean folderMade = chromDir.mkdir();
		
		if (!chromDir.exists() && !folderMade) {
			spacer = "_";
		}
		
		//-- Read input files and split into chromosome files--
		if (!sc.hasNextLine()) {
			sc.close();
			return null;
		}
		
		String line = null;
		String chrom = null;
		while (chrom == null && sc.hasNextLine()) {
			line = sc.nextLine();
			chrom = parseLineChrom(line);
		}
		
		File cFile = new File(chromDirName + spacer + chrom + ".bed");
		chromFiles.add(cFile);
		PrintWriter pw = new PrintWriter(cFile);
		pw.println(line);
		
		//separate files
		while(sc.hasNextLine()) { 
			String nChrom = null;
			while (nChrom == null && sc.hasNextLine()) {
				line = sc.nextLine();
				nChrom = parseLineChrom(line);
			}
			if (!chrom.equals(nChrom)) {
				pw.close();
				chrom = nChrom;
				cFile = new File(chromDirName + spacer + nChrom + ".bed");
				chromFiles.add(cFile);
				pw = new PrintWriter(cFile);
			}
			pw.println(line);
		}
		
		pw.close();
		sc.close();
		
		return chromFiles;
	}
	
	/**
	 * Checks if the given input line is able to be parsed
	 * @param line - line to be parsed
	 * @return the chromosome "chr_" given by the line, null if unable to be fully parsed 
	 */
	public static String parseLineChrom(String line) {
		String[] ra = line.split("\\s+");

		if (ra[0].startsWith("chr")) {
			return ra[0];
		}
		
		try {
			Integer.parseInt(ra[1]);
			Integer.parseInt(ra[2]);
			Integer.parseInt(ra[3]);
		} catch (NumberFormatException e) {
			return null;
		} catch (IndexOutOfBoundsException e) {
			return null;
		}
		return null;
	}
	
	/**
	 * Concatenate output chromosome files into one genome file (in order that it was given in)
	 * @param outChrFiles - ArrayList of output chromosome files
	 * @throws FileNotFoundException - if unable to find a chromosome file, or if unable to write to file
	 */
	public static void mergeWholeChrFile(ArrayList<File> outChrFiles) throws FileNotFoundException {
		String wChrFileName = chromDirName.substring(0, chromDirName.lastIndexOf("_ChromFiles"))
				+ "_out" + spacer + outBaseName + "_allChr.bed";
		File wChrFile = new File (wChrFileName);
		PrintWriter wholePW = null;
		try {
			wholePW = new PrintWriter(wChrFile);
		} catch (FileNotFoundException e) {
			throw new FileNotFoundException("File not found: " + wChrFileName);
		}
		
		/*wholePW.println("track type=bedGraph name=\"" + wChrFileName + "\"" + " description=\"" + wChrFileName + "\" "
				+ "visibility=full autoScale=Off alwaysZero=On maxHeightPixels=128:30:11 viewLimits=0:1"); //header */
		
		for(File cFile : outChrFiles) {
			Scanner sc = new Scanner(cFile);
			/*sc.nextLine(); //skip chromosome file headers*/
			while(sc.hasNextLine()) {
				wholePW.println(sc.nextLine());
			}
			sc.close();
		}
		
		wholePW.close();
		System.out.println("Merging completed");
	}
	
}
