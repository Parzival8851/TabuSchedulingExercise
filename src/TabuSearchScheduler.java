import java.util.*;

/**
 * TabuSearchScheduler implements a Tabu Search metaheuristic to solve a scheduling problem
 * where the objective is to minimize the total weighted tardiness of jobs.
 *
 * <p>This implementation uses:
 * <ul>
 *   <li>Adaptive tabu tenure and jump size for neighborhood exploration based on adjacent and non-adjacent swaps</li>
 *   <li>A dynamic aspiration threshold that relaxes when improvements are scarce</li>
 *   <li>A Long-Term Memory (LTM) frequency-based penalty to discourage overused moves</li>
 * </ul>
 *
 * <p>All parameters influencing the Tabu Search behavior are defined as final constants
 * at the beginning of the program, grouped by role.
 *
 * @author Matteo Boniotti - 731246
 */
public class TabuSearchScheduler {
	// ============================ Tabu Search Parameters ============================

	// ---------------------------- General Parameters ----------------------------
	/** Maximum factorial value used to cap the maximum number of iterations */
	public static final int MAX_FACTORIAL = 10_000;

	// ---------------------------- Tenure Parameters ----------------------------
	/** Minimum tabu tenure fraction (relative to number of jobs) */
	public static final double MIN_TENURE_FRAC = 0.05;
	/** Initial tabu tenure fraction (relative to number of jobs) */
	public static final double INITIAL_TENURE_FRAC = 0.15;
	/** Maximum tabu tenure fraction (relative to number of jobs) */
	public static final double MAX_TENURE_FRAC = 0.75;
	/** Fraction used to calculate the tenure decay threshold relative to maximum iterations */
	public static final double TENURE_DECAY_THRESHOLD = 0.5;

	// ---------------------------- No Improvement Parameters ----------------------------
	/** Fraction used to calculate the no-improvement threshold relative to maximum iterations */
	public static final int NO_IMPROVEMENT_THRESHOLD_FRAC = 100;
	/** Minimum number of iterations without improvement before adapting parameters */
	public static final int MIN_NO_IMPROVEMENT_THRESHOLD = 50;

	// ---------------------------- Aspiration Criterion Parameters ----------------------------
	/** Fraction used to calculate the aspiration threshold relative to best solution value */
	public static final double BASE_ASPIRATION = 0.05;
	/** Minimum aspiration threshold */
	public static final double MIN_ASPIRATION = 0.01;
	/** Fraction used to adapt the aspiration threshold */
	public static final double FRAC_ASPIRATION_ADAPT = 0.04;

	// ---------------------------- Long-Term Memory (LTM) Parameters ----------------------------
	/** Coefficient used for penalizing move frequency in the LTM frequency-based mechanism */
	public static final double PENALTY_COEFFICIENT = 0.1;

	// ---------------------------- Quick Evaluation Parameters ----------------------------
	/** Fraction used to calculate the quick evaluation threshold relative to maximum iterations */
	public static final double FRAC_QUICK_EVAL = 0.5;

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
		int maxIterations = factorial(size, MAX_FACTORIAL);

		// Tabu tenure settings (absolute values)
		// Tabu tenure is defined as 1 <= min <= initial <= max, relative to the number of jobs.
		// The initial tenure is set between min and max, based on the number of jobs, to allow for an exploration phase at the beginning
		// The tabu tenure is not already set to the min or the max value to allow for a middle phase where the algorithm can adapt to the problem.
		int maxTabuTenure = (int) (jobs.size() * MAX_TENURE_FRAC);
		int minTabuTenure = (int) Math.max(size * MIN_TENURE_FRAC, 1);
		int initialTabuTenure = (int) Math.max(minTabuTenure, Math.min(size * INITIAL_TENURE_FRAC, maxTabuTenure));

		// Calculate the no-improvement threshold: the number of iterations without improvement before adaptation (minimum threshold enforced)
		// The threshold enables the algorithm to adapt its parameters when stuck in a local minimum, by increasing the tabu tenure and jump size.
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
		// INITIAL SOLUTION: Use the Earliest Due Date (EDD) rule - sort jobs by due date.
		Solution currentSolution = new Solution(jobs);
		currentSolution.sortByDueDate();

		System.out.println("EDD Initial Solution:");
		System.out.println(currentSolution);

		// Global variable to store the best solution found.
		Solution bestSolution = new Solution(currentSolution);
		double bestObjective = bestSolution.getObjective();

		// Iteration tracking.
		int bestIteration = 0;
		int iteration = 0;
		int currentTabuTenure = initialTabuTenure;
		int nonImprovementCount = 0;

		// Adaptive jump size for non-adjacent swaps.
		int minJumpSize = 1;
		int maxJumpSize = currentSolution.size() - 1;
		int currentJumpSize = minJumpSize;

		// Tabu list: move string -> iteration when it expires.
		Map<String, Integer> tabuList = new HashMap<>();
		// LTM frequency map: move string -> frequency count.
		Map<String, Integer> frequencyMap = new HashMap<>();

		// Main Tabu Search loop.
		while (iteration < maxIterations) {
			iteration++;
			List<Candidate> allCandidateList = new ArrayList<>(); // Store all generated candidates.
			List<Candidate> validCandidateList = new ArrayList<>(); // Store valid candidates (not tabu or meeting aspiration criterion).

			// --------------------- Dynamic Aspiration Threshold ---------------------
			// Base threshold is 5%; relax it as non-improvement increases even so slight.
			double dynamicAspiration = BASE_ASPIRATION - ((double) nonImprovementCount / noImprovementThreshold) * FRAC_ASPIRATION_ADAPT;
			dynamicAspiration = Math.max(dynamicAspiration, MIN_ASPIRATION);

			// --------------------- Candidate Generation ---------------------
			/*
			 * Defines the neighborhood N(x, jump) by swapping component i with component j = i + k * jump,
			 * for k = 1, 2, ..., as long as j <= n.
			 *
			 * Given a vector x = (x1, x2, ..., xn) and a positive integer jump ∈ N^+,
			 * the neighborhood is defined as:
			 *
			 * N(x, jump) = { y ∈ R^n : ∃ i ∈ {1, ..., n}, k ∈ N^+ such that
			 * j = i + k * jump ≤ n and
			 * y = (x1, ..., x_{i-1}, x_j, x_{i+1}, ..., x_{j-1}, x_i, x_{j+1}, ..., x_n) }
			 *
			 * Here:
			 * - i: the index of the component to be swapped.
			 * - k: the step count, determining how far to jump.
			 * - j = i + k * jump: the index of the component with which x_i is swapped.
			 *
			 * This neighborhood explores permutations obtained by periodic jumps
			 * along the vector based on the given jump size.
			 */
			for (int i = 0; i < currentSolution.size() - currentJumpSize; i++) {
				for (int j = i + currentJumpSize; j < currentSolution.size(); j += currentJumpSize) {
					// Create a neighbor solution by swapping jobs at positions i and j.
					Solution neighbor = new Solution(currentSolution);
					neighbor.swap(i, j);
					// Evaluate the neighbor solution and compute its quick evaluation (maximum tardiness).
					double neighborObjective = neighbor.getObjective();
					double quickEval = neighbor.getMaxTardiness();

					// Build a normalized move key (e.g., "SWAP:1:3")
					// The move is always represented as the smaller ID first because the jump are always forward.
					int id1 = Math.min(currentSolution.getJobIdAt(i), currentSolution.getJobIdAt(j));
					int id2 = Math.max(currentSolution.getJobIdAt(i), currentSolution.getJobIdAt(j));
					String moveKey = "SWAP:" + id1 + ":" + id2;

					// Create a candidate solution as the initial solution with the proposed move.
					Candidate candidate = new Candidate(neighbor, moveKey, neighborObjective, quickEval);
					allCandidateList.add(candidate);

					// Check tabu status and incorporate LTM frequency penalty.
					boolean isTabu = tabuList.containsKey(moveKey) && tabuList.get(moveKey) > iteration;
					int freq = frequencyMap.getOrDefault(moveKey, 0);
					// LTM here is used as a penalization for already used moves so the search prefers new moves in its exploration
					double adjustedObj = neighborObjective + PENALTY_COEFFICIENT * freq;

					// Apply dynamic aspiration: allow candidate if not tabu, or if it meets the aspiration threshold.
					// Aspiration threshold is a percentage of the best objective found so far: if the candidate is better than
					// the best objective by this percentage, it is accepted even if it is tabu.
					if (!isTabu || adjustedObj <= bestObjective * (1 - dynamicAspiration)) {
						validCandidateList.add(candidate);
					}
				}
			}

			// Fallback: if no candidate qualifies under tabu/aspiration, use all generated candidates.
			// This ensures that the algorithm always makes a move, even if all candidates are tabu or not aspiring.
			// The tabu list is ignored in this case, and a soft reset is applied by removing entries older than a certain threshold.
			List<Candidate> candidateList;

			int finalIteration = iteration;
			if (validCandidateList.isEmpty()) {
				// Calculate the decay threshold dynamically based on the current tabu tenure
				int decayThreshold = Math.max(1, (int) (currentTabuTenure * TENURE_DECAY_THRESHOLD));

				// Remove tabu entries older than the calculated threshold
				tabuList.entrySet().removeIf(entry -> entry.getValue() < finalIteration + decayThreshold);

				// Use all candidates as fallback
				candidateList = allCandidateList;
			} else {
				candidateList = validCandidateList;
			}

			// --------------------- Two-Phase Candidate Selection ---------------------
			// Phase 1: Sort by quick evaluation (maximum tardiness).
			candidateList.sort(Comparator.comparingDouble(c -> c.quickEvaluation));
			// Select a fraction of the candidates for full evaluation.
			int subsetSize = (int) Math.max(1, candidateList.size() * FRAC_QUICK_EVAL);
			List<Candidate> candidateSubset = candidateList.subList(0, subsetSize);
			// Phase 2: Select the candidate with the lowest full objective.
			// In this second phase we could add some randomness with an RCL list to avoid getting stuck in local optima.
			// In this implementation we simply select the best candidate.
			Candidate bestCandidate = Collections.min(candidateSubset, Comparator.comparingDouble(c -> c.objectiveValue));

			// Update the current solution.
			currentSolution = bestCandidate.solution;
			double currentObjective = bestCandidate.objectiveValue;

			// --------------------- Update Long-Term Memory (LTM) ---------------------
			// Increment the frequency count for the performed move to reduce its likelihood of being used again
			// promoting better exploration of the solution space.
			frequencyMap.put(bestCandidate.moveKey, frequencyMap.getOrDefault(bestCandidate.moveKey, 0) + 1);


			// --------------------- Adaptive Parameter Updates ---------------------
			if (currentObjective < bestObjective) {
				// Improvement: update the best solution and reset parameters to intensify the search.
				bestSolution = new Solution(currentSolution);
				bestObjective = currentObjective;
				bestIteration = iteration;
				nonImprovementCount = 0;
				currentTabuTenure = minTabuTenure;
				currentJumpSize = minJumpSize;
			} else {
				// No improvement: increase non-improvement count and possibly adapt parameters to diversify the search.
				// This is a simple strategy to adapt the tabu tenure and jump size based on the number of non-improving iterations.
				nonImprovementCount++;
				if (nonImprovementCount >= noImprovementThreshold) {
					currentTabuTenure = Math.min(currentTabuTenure + 1, maxTabuTenure);
					currentJumpSize = Math.min(currentJumpSize + 1, maxJumpSize);
					nonImprovementCount = 0;
				}
			}

			// --------------------- Update Tabu List ---------------------
			// Mark the move as tabu for the next 'currentTabuTenure' iterations.
			// The tabu list stores direct moves, but it can be easily modified to include reverse moves
			// by adjusting the indexes in the moveKey.
			// To track backward swaps, a new moveKey format like "REVERSE:id1:id2" could be introduced.
			tabuList.put(bestCandidate.moveKey, iteration + currentTabuTenure);
			tabuList.entrySet().removeIf(entry -> entry.getValue() <= finalIteration);
		}

		// --------------------- Return the Best Solution ---------------------
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
	public static int factorial(int n, int maxValue) {
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
 * <p>
 * This class caches the computed total weighted tardiness and maximum tardiness.
 * It uses a dirty flag to recompute these values only when the schedule has been modified.
 */
class Solution {
	List<Job> schedule;

	// Cached objective values.
	private double cachedObjective;
	private double cachedMaxTardiness;

	// Flag to indicate if the schedule has been modified and the cached values are no longer valid.
	private boolean isDirty;

	/**
	 * Constructs a Solution using the given list of jobs.
	 *
	 * @param jobs the list of jobs.
	 */
	public Solution(List<Job> jobs) {
		this.schedule = new ArrayList<>(jobs);
		this.isDirty = true; // Initially, we need to compute the objectives.
	}

	/**
	 * Copy constructor.
	 *
	 * @param other the Solution to copy.
	 */
	public Solution(Solution other) {
		this.schedule = new ArrayList<>(other.schedule);
		this.isDirty = true; // New copy should compute objectives anew.
	}

	/**
	 * Sorts the schedule by due date (Earliest Due Date first) and marks the solution as modified.
	 */
	public void sortByDueDate() {
		schedule.sort(Comparator.comparingInt(job -> job.dueDate));
		isDirty = true;
	}

	/**
	 * Returns the total weighted tardiness of the schedule.
	 * Recomputes the value only if the schedule was modified.
	 *
	 * @return the total weighted tardiness.
	 */
	public double getObjective() {
		if (isDirty) {
			computeObjectives();
		}
		return cachedObjective;
	}

	/**
	 * Returns the maximum tardiness among all jobs in the schedule.
	 * Recomputes the value only if the schedule was modified.
	 *
	 * @return the maximum tardiness.
	 */
	public double getMaxTardiness() {
		if (isDirty) {
			computeObjectives();
		}
		return cachedMaxTardiness;
	}

	/**
	 * Computes both the total weighted tardiness and maximum tardiness in one pass.
	 * Updates the cached values and resets the dirty flag.
	 */
	public void computeObjectives() {
		double total = 0;
		double maxTardiness = 0;
		int currentTime = 0;
		for (Job job : schedule) {
			currentTime += job.processingTime;
			double tardiness = Math.max(0, currentTime - job.dueDate);
			total += job.weight * tardiness;
			if (tardiness > maxTardiness) {
				maxTardiness = tardiness;
			}
		}
		cachedObjective = total;
		cachedMaxTardiness = maxTardiness;
		isDirty = false;
	}

	/**
	 * Swaps the jobs at indices i and j and marks the solution as modified.
	 *
	 * @param i the first index.
	 * @param j the second index.
	 */
	public void swap(int i, int j) {
		Collections.swap(schedule, i, j);
		isDirty = true;
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
	 * Returns the job ID at the specified index.
	 *
	 * @param index the index.
	 * @return the job ID.
	 */
	public int getJobIdAt(int index) {
		return schedule.get(index).id;
	}

	@Override
	public String toString() {
		// Ensure the objectives are up-to-date.
		if (isDirty) {
			computeObjectives();
		}
		StringBuilder sb = new StringBuilder();
		sb.append("Schedule: ");
		for (Job job : schedule) {
			sb.append(job).append(" ");
		}
		sb.append("\nTotal Weighted Tardiness: ").append(cachedObjective);
		return sb.toString();
	}
}

