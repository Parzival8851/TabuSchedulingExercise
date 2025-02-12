import java.util.*;

/**
 * TabuSearchScheduler implements a Tabu Search metaheuristic to solve a scheduling problem
 * where the objective is to minimize the total weighted tardiness of jobs.
 *
 * <p>This implementation uses:
 * <ul>
 *   <li>Adaptive tabu tenure and jump size for neighborhood exploration</li>
 *   <li>A dynamic aspiration threshold that relaxes when improvements are scarce</li>
 *   <li>A Long-Term Memory (LTM) frequency-based penalty to discourage overused moves</li>
 *   <li>A dedicated {@link Solution} class to encapsulate the schedule operations</li>
 * </ul>
 *
 * <p>All parameters influencing the Tabu Search behavior are defined as final constants
 * at the beginning of the program, grouped by role.
 *
 * @author
 */
public class TabuSearchScheduler {

	// ============================ General Parameters ============================
	/** Maximum factorial value used to cap the maximum number of iterations */
	public static final int MAX_FACTORIAL = 10_000;

	// ============================ Tabu Search Parameters ============================
	/** Minimum tabu tenure fraction (relative to number of jobs) */
	public static final double MIN_TENURE_FRAC = 0.005;
	/** Initial tabu tenure fraction (relative to number of jobs) */
	public static final double INITIAL_TENURE_FRAC = 0.15;
	/** Maximum tabu tenure fraction (relative to number of jobs) */
	public static final double MAX_TENURE_FRAC = 0.75;

	/** Fraction used to calculate the no-improvement threshold relative to maximum iterations */
	public static final int NO_IMPROVEMENT_THRESHOLD_FRAC = 100;
	/** Minimum number of iterations without improvement before adapting parameters */
	public static final int MIN_NO_IMPROVEMENT_THRESHOLD = 50;

	// ============================ Long-Term Memory (LTM) Parameters ============================
	/** Coefficient used for penalizing move frequency in the LTM frequency-based mechanism */
	public static final double PENALTY_COEFFICIENT = 0.1;

	/**
	 * Main method.
	 *
	 * @param args Command line arguments (not used)
	 */
	public static void main(String[] args) {
		int w = 1;
		// Define jobs with processing time (p), due date (d), and weight (w)
		List<Job> jobs = Arrays.asList(
				new Job(1, 6, 9, w),
				new Job(2, 4, 12, w),
				new Job(3, 8, 15, w),
				new Job(4, 2, 8, w),
				new Job(5, 10, 20, w),
				new Job(6, 3, 22, w)
		);

		// --------------------- Parameter Setup ---------------------
		int size = jobs.size();
		// Maximum iterations is capped by the factorial of the number of jobs or MAX_FACTORIAL.
		int maxIterations = fattoriale(size, MAX_FACTORIAL);

		// Tabu tenure settings (absolute values)
		int maxTabuTenure = (int) (jobs.size() * MAX_TENURE_FRAC);
		int minTabuTenure = (int) Math.max(size * MIN_TENURE_FRAC, 1);
		int initialTabuTenure = (int) Math.max(minTabuTenure, Math.min(size * INITIAL_TENURE_FRAC, maxTabuTenure));

		// Calculate the no-improvement threshold: the number of iterations without improvement before adaptation.
		int noImprovementThreshold = Math.max(maxIterations / NO_IMPROVEMENT_THRESHOLD_FRAC, MIN_NO_IMPROVEMENT_THRESHOLD);

		// --------------------- Run Tabu Search ---------------------
		Solution bestSolution = tabuSearch(jobs, maxIterations, initialTabuTenure, minTabuTenure, maxTabuTenure, noImprovementThreshold);

		System.out.println("Best solution found:");
		System.out.println(bestSolution);
	}

	/**
	 * Performs the Tabu Search metaheuristic to minimize total weighted tardiness.
	 *
	 * <p>This method generates candidate moves using non-adjacent swaps and uses adaptive parameters,
	 * a dynamic aspiration threshold, and a Long-Term Memory (LTM) mechanism.
	 *
	 * @param jobs the list of jobs to schedule.
	 * @param maxIterations the maximum number of iterations.
	 * @param initialTabuTenure the initial tabu tenure.
	 * @param minTabuTenure the minimum tabu tenure.
	 * @param maxTabuTenure the maximum tabu tenure.
	 * @param noImprovementThreshold the threshold (in iterations) without improvement before adapting parameters.
	 * @return the best solution found.
	 */
	public static Solution tabuSearch(List<Job> jobs, int maxIterations,
									  int initialTabuTenure, int minTabuTenure, int maxTabuTenure, int noImprovementThreshold) {
		// INITIAL SOLUTION: Use the Earliest Due Date (EDD) rule.
		Solution currentSolution = new Solution(jobs);
		currentSolution.sortByDueDate();

		Solution bestSolution = new Solution(currentSolution);
		double bestObjective = bestSolution.getObjective();

		System.out.println("EDD Initial Solution:");
		System.out.println(currentSolution);

		// Iteration tracking.
		int bestIteration = 0;
		int iteration = 0;
		int currentTabuTenure = initialTabuTenure;
		int nonImprovementCount = 0;

		// Adaptive jump size for non-adjacent swaps.
		int minJumpSize = 1;
		int maxJumpSize = currentSolution.size() - 1;
		int currentJumpSize = minJumpSize;

		// Tabu list: maps move keys (String) to iteration numbers until which the move is tabu.
		Map<String, Integer> tabuList = new HashMap<>();
		// LTM frequency map: counts the frequency of each move.
		Map<String, Integer> frequencyMap = new HashMap<>();

		// Main Tabu Search loop.
		while (iteration < maxIterations) {
			iteration++;
			List<Candidate> allCandidateList = new ArrayList<>();
			List<Candidate> validCandidateList = new ArrayList<>();

			// --------------------- Dynamic Aspiration Threshold ---------------------
			// Base threshold is 5%; relax it as non-improvement increases.
			double baseAspiration = 0.05;
			double dynamicAspiration = baseAspiration - ((double) nonImprovementCount / noImprovementThreshold) * 0.04;
			dynamicAspiration = Math.max(dynamicAspiration, 0.01);

			// --------------------- Candidate Generation ---------------------
			// Generate candidate moves using non-adjacent swaps.
			for (int i = 0; i < currentSolution.size() - currentJumpSize; i++) {
				// For each i, generate moves for j = i+currentJumpSize, i+2*currentJumpSize, etc.
				for (int j = i + currentJumpSize; j < currentSolution.size(); j += currentJumpSize) {
					// Create a neighbor solution by swapping jobs at positions i and j.
					Solution neighbor = new Solution(currentSolution);
					neighbor.swap(i, j);
					double neighborObjective = neighbor.getObjective();
					double quickEval = neighbor.getMaxTardiness();

					// Build a normalized move key (e.g., "SWAP:1:3").
					int id1 = Math.min(currentSolution.getJobIdAt(i), currentSolution.getJobIdAt(j));
					int id2 = Math.max(currentSolution.getJobIdAt(i), currentSolution.getJobIdAt(j));
					String moveKey = "SWAP:" + id1 + ":" + id2;

					Candidate candidate = new Candidate(neighbor, moveKey, neighborObjective, quickEval);
					allCandidateList.add(candidate);

					// Check tabu status and incorporate LTM frequency penalty.
					boolean isTabu = tabuList.containsKey(moveKey) && tabuList.get(moveKey) > iteration;
					int freq = frequencyMap.getOrDefault(moveKey, 0);
					double adjustedObj = neighborObjective + PENALTY_COEFFICIENT * freq;

					// Apply dynamic aspiration: allow candidate if not tabu, or if it meets the aspiration threshold.
					if (!isTabu || adjustedObj <= bestObjective * (1 - dynamicAspiration)) {
						validCandidateList.add(candidate);
					}
				}
			}

			// Fallback: if no candidate qualifies under tabu/aspiration, use all generated candidates.
			List<Candidate> candidateList = validCandidateList.isEmpty() ? allCandidateList : validCandidateList;

			// --------------------- Two-Phase Candidate Selection ---------------------
			// Phase 1: Sort by quick evaluation (maximum tardiness).
			candidateList.sort(Comparator.comparingDouble(c -> c.quickEvaluation));
			// Select the best 50% (by quick evaluation).
			int subsetSize = Math.max(1, candidateList.size() / 2);
			List<Candidate> candidateSubset = candidateList.subList(0, subsetSize);
			// Phase 2: Select the candidate with the lowest full objective.
			Candidate bestCandidate = Collections.min(candidateSubset, Comparator.comparingDouble(c -> c.objectiveValue));

			// Update the current solution.
			currentSolution = bestCandidate.solution;
			double currentObjective = bestCandidate.objectiveValue;

			// --------------------- Update Long-Term Memory (LTM) ---------------------
			// Increment the frequency count for the move performed.
			frequencyMap.put(bestCandidate.moveKey, frequencyMap.getOrDefault(bestCandidate.moveKey, 0) + 1);

			// --------------------- Adaptive Parameter Updates ---------------------
			if (currentObjective < bestObjective) {
				// Improvement: update best solution and reset parameters.
				bestSolution = new Solution(currentSolution);
				bestObjective = currentObjective;
				bestIteration = iteration;
				nonImprovementCount = 0;
				currentTabuTenure = minTabuTenure;
				currentJumpSize = minJumpSize; // Reset jump size on improvement.
			} else {
				// No improvement: increase non-improvement count and possibly adapt parameters.
				nonImprovementCount++;
				if (nonImprovementCount >= noImprovementThreshold) {
					currentTabuTenure = Math.min(currentTabuTenure + 1, maxTabuTenure);
					currentJumpSize = Math.min(currentJumpSize + 1, maxJumpSize);
					nonImprovementCount = 0;
				}
			}

			// --------------------- Update Tabu List ---------------------
			// Mark the move as tabu for the next 'currentTabuTenure' iterations.
			tabuList.put(bestCandidate.moveKey, iteration + currentTabuTenure);
			int finalIteration = iteration;
			tabuList.entrySet().removeIf(entry -> entry.getValue() <= finalIteration);
		}

		System.out.println("Best solution found at iteration: " + bestIteration);
		return bestSolution;
	}

	/**
	 * Computes the factorial of n, capped by maxValue to prevent excessive iterations.
	 *
	 * @param n the input number.
	 * @param maxValue the maximum allowed value.
	 * @return the factorial of n or the current result if it exceeds maxValue.
	 */
	public static int fattoriale(int n, int maxValue) {
		if (n <= 1) {
			return 1;
		}
		int result = n;
		for (int i = n - 1; i >= 2; i--) {
			if (result > maxValue / i) { // Prevent fictitious overflow.
				return result;
			}
			result *= i;
		}
		return result;
	}

	/**
	 * Computes the total weighted tardiness for a sequence of jobs.
	 *
	 * @param sequence the list of jobs.
	 * @return the total weighted tardiness.
	 */
	public static double calculateTotalWeightedTardiness(List<Job> sequence) {
		double totalWeightedTardiness = 0;
		double currentTime = 0;
		for (Job job : sequence) {
			currentTime += job.processingTime;
			double tardiness = Math.max(0, currentTime - job.dueDate);
			totalWeightedTardiness += job.weight * tardiness;
		}
		return totalWeightedTardiness;
	}

	/**
	 * Computes the maximum tardiness in a sequence of jobs.
	 *
	 * @param sequence the list of jobs.
	 * @return the maximum tardiness.
	 */
	public static double calculateMaxTardiness(List<Job> sequence) {
		double maxTardiness = 0;
		double currentTime = 0;
		for (Job job : sequence) {
			currentTime += job.processingTime;
			double tardiness = Math.max(0, currentTime - job.dueDate);
			if (tardiness > maxTardiness) {
				maxTardiness = tardiness;
			}
		}
		return maxTardiness;
	}
}

// ============================ Supporting Classes ============================

/**
 * Represents a job with an ID, processing time, due date, and weight (penalty for tardiness).
 */
class Job {
	int id;
	int processingTime;
	int dueDate;
	int weight;

	/**
	 * Constructs a Job instance.
	 *
	 * @param id the job ID.
	 * @param processingTime the processing time.
	 * @param dueDate the due date.
	 * @param weight the weight (penalty multiplier).
	 */
	public Job(int id, int processingTime, int dueDate, int weight) {
		this.id = id;
		this.processingTime = processingTime;
		this.dueDate = dueDate;
		this.weight = weight;
	}

	@Override
	public String toString() {
		return "Job" + id;
	}
}

/**
 * Represents a move between two jobs in a schedule.
 * The move is defined by two job IDs, stored in sorted order.
 */
class Move {
	int jobId1;
	int jobId2;

	/**
	 * Constructs a Move between two job IDs.
	 *
	 * @param jobId1 the first job ID.
	 * @param jobId2 the second job ID.
	 */
	public Move(int jobId1, int jobId2) {
		// Normalize the order.
		this.jobId1 = Math.min(jobId1, jobId2);
		this.jobId2 = Math.max(jobId1, jobId2);
	}

	@Override
	public int hashCode() {
		return Objects.hash(jobId1, jobId2);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj instanceof Move) {
			Move other = (Move) obj;
			return this.jobId1 == other.jobId1 && this.jobId2 == other.jobId2;
		}
		return false;
	}
}

/**
 * Holds a candidate solution along with the move key used to generate it,
 * its full objective value (total weighted tardiness), and a quick evaluation (maximum tardiness).
 */
class Candidate {
	Solution solution;
	String moveKey;
	double objectiveValue;    // Full objective value
	double quickEvaluation;   // Quick evaluation value

	/**
	 * Constructs a Candidate.
	 *
	 * @param solution the candidate solution.
	 * @param moveKey the key representing the move (e.g., "SWAP:1:3").
	 * @param objectiveValue the full objective value.
	 * @param quickEvaluation the quick evaluation (maximum tardiness).
	 */
	public Candidate(Solution solution, String moveKey, double objectiveValue, double quickEvaluation) {
		this.solution = solution;
		this.moveKey = moveKey;
		this.objectiveValue = objectiveValue;
		this.quickEvaluation = quickEvaluation;
	}
}

/**
 * Encapsulates a schedule (list of jobs) and provides operations for evaluation and manipulation,
 * including sorting by due date, swapping jobs, and computing objective values.
 */
class Solution {
	List<Job> schedule;

	/**
	 * Constructs a Solution using the given list of jobs.
	 *
	 * @param jobs the list of jobs.
	 */
	public Solution(List<Job> jobs) {
		this.schedule = new ArrayList<>(jobs);
	}

	/**
	 * Copy constructor.
	 *
	 * @param other the Solution to copy.
	 */
	public Solution(Solution other) {
		this.schedule = new ArrayList<>(other.schedule);
	}

	/**
	 * Sorts the schedule by due date (Earliest Due Date first).
	 */
	public void sortByDueDate() {
		schedule.sort(Comparator.comparingInt(job -> job.dueDate));
	}

	/**
	 * Computes the total weighted tardiness of the schedule.
	 *
	 * @return the total weighted tardiness.
	 */
	public double getObjective() {
		double total = 0;
		int currentTime = 0;
		for (Job job : schedule) {
			currentTime += job.processingTime;
			total += job.weight * Math.max(0, currentTime - job.dueDate);
		}
		return total;
	}

	/**
	 * Computes the maximum tardiness among all jobs in the schedule.
	 *
	 * @return the maximum tardiness.
	 */
	public double getMaxTardiness() {
		double maxTardiness = 0;
		int currentTime = 0;
		for (Job job : schedule) {
			currentTime += job.processingTime;
			double tardiness = Math.max(0, currentTime - job.dueDate);
			if (tardiness > maxTardiness) {
				maxTardiness = tardiness;
			}
		}
		return maxTardiness;
	}

	/**
	 * Computes the tardiness for the job at the specified index.
	 *
	 * @param index the index of the job.
	 * @return the tardiness of the job at that index.
	 */
	public double getJobTardiness(int index) {
		int currentTime = 0;
		for (int i = 0; i <= index; i++) {
			currentTime += schedule.get(i).processingTime;
		}
		Job job = schedule.get(index);
		return Math.max(0, currentTime - job.dueDate);
	}

	/**
	 * Returns the number of jobs in the schedule.
	 *
	 * @return the size of the schedule.
	 */
	public int size() {
		return schedule.size();
	}

	/**
	 * Swaps the jobs at indices i and j.
	 *
	 * @param i the first index.
	 * @param j the second index.
	 */
	public void swap(int i, int j) {
		Collections.swap(schedule, i, j);
	}

	/**
	 * Returns the job ID at the specified index.
	 *
	 * @param index the index.
	 * @return the job ID.
	 */
	public int getJobIdAt(int index) {
		return schedule.get(index).id;
	}

	/**
	 * Removes and returns the job at the specified index.
	 *
	 * @param index the index.
	 * @return the removed job.
	 */
	public Job removeJob(int index) {
		return schedule.remove(index);
	}

	/**
	 * Inserts a job at the specified index.
	 *
	 * @param index the index.
	 * @param job the job to insert.
	 */
	public void insertJob(int index, Job job) {
		schedule.add(index, job);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Schedule: ");
		for (Job job : schedule) {
			sb.append(job).append(" ");
		}
		sb.append("\nTotal Weighted Tardiness: ").append(getObjective());
		return sb.toString();
	}
}
