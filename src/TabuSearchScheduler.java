import java.util.*;

public class TabuSearchScheduler {

	// --- Parameters ---
	public static final int MAX_FACTORIAL = 10_000;

	public static final double MIN_TENURE_FRAC = 0.005;
	public static final double INITIAL_TENURE_FRAC = 0.15;
	public static final double MAX_TENURE_FRAC = 0.75;

	public static final int NO_IMPROVEMENT_THRESHOLD_FRAC = 100;
	public static final int MIN_NO_IMPROVEMENT_THRESHOLD = 50;

	// Coefficient for long-term memory penalty
	public static final double PENALTY_COEFFICIENT = 0.1;

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

		// Setup parameters.
		int size = jobs.size();
		int maxIterations = fattoriale(size, MAX_FACTORIAL); // maxIterations = min(n!, MAX_FACTORIAL)

		int maxTabuTenure = (int)(jobs.size() * MAX_TENURE_FRAC);
		int minTabuTenure = (int)Math.max(size * MIN_TENURE_FRAC, 1);
		int initialTabuTenure = (int)Math.max(minTabuTenure, Math.min(size * INITIAL_TENURE_FRAC, maxTabuTenure));
		int noImprovementThreshold = Math.max(maxIterations / NO_IMPROVEMENT_THRESHOLD_FRAC, MIN_NO_IMPROVEMENT_THRESHOLD);

		// Run the enhanced Tabu Search.
		Solution bestSolution = tabuSearch(jobs, maxIterations,
				initialTabuTenure, minTabuTenure, maxTabuTenure, noImprovementThreshold);

		System.out.println("Best solution found:");
		System.out.println(bestSolution);
	}

	public static Solution tabuSearch(List<Job> jobs, int maxIterations,
									  int initialTabuTenure, int minTabuTenure, int maxTabuTenure, int noImprovementThreshold) {

		// INITIAL SOLUTION: use EDD (Earliest Due Date) rule.
		Solution currentSolution = new Solution(jobs);
		currentSolution.sortByDueDate();

		Solution bestSolution = new Solution(currentSolution);
		double bestObjective = bestSolution.getObjective();
		int bestIteration = 0;

		// Tabu list: keys are strings that encode the move (e.g. "SWAP:1:2" or "INSERT:3:5->2").
		Map<String, Integer> tabuList = new HashMap<>();
		// Long Term Memory: frequency map to count how many times a move is performed.
		Map<String, Integer> frequencyMap = new HashMap<>();

		int iteration = 0;
		int currentTabuTenure = initialTabuTenure;
		int nonImprovementCount = 0;

		// Adaptive jump size parameters.
		int minJumpSize = 1;
		int maxJumpSize = currentSolution.size() - 1;
		int currentJumpSize = minJumpSize;

		while (iteration < maxIterations) {
			iteration++;
			List<Candidate> allCandidates = new ArrayList<>();

			// Compute average tardiness of the current solution (used for “jumping back” criteria).
			double totalTardiness = 0;
			for (int i = 0; i < currentSolution.size(); i++) {
				totalTardiness += currentSolution.getJobTardiness(i);
			}
			double avgTardiness = totalTardiness / currentSolution.size();

			// --- Generate forward swap candidates ---
			// For each valid index i, swap with job at index i+currentJumpSize.
			for (int i = 0; i < currentSolution.size() - currentJumpSize; i++) {
				int j = i + currentJumpSize;
				Solution neighbor = new Solution(currentSolution);
				neighbor.swap(i, j);
				double obj = neighbor.getObjective();
				double quickEval = neighbor.getMaxTardiness();
				// Create a move key for the swap (order is normalized).
				String moveKey = "SWAP:"
						+ Math.min(neighbor.getJobIdAt(i), neighbor.getJobIdAt(j)) + ":"
						+ Math.max(neighbor.getJobIdAt(i), neighbor.getJobIdAt(j));
				Candidate cand = new Candidate(neighbor, moveKey, obj, quickEval, false);
				allCandidates.add(cand);
			}

			// --- Generate "jumping back" (insertion) candidates ---
			// For jobs that are in a "high" index and have tardiness above average, try moving them earlier.
			for (int i = currentJumpSize; i < currentSolution.size(); i++) {
				double tardiness = currentSolution.getJobTardiness(i);
				if (tardiness > avgTardiness) {  // job with "high" tardiness
					int newPos = i - currentJumpSize;
					if (newPos >= 0) {
						Solution neighbor = new Solution(currentSolution);
						Job job = neighbor.removeJob(i);
						neighbor.insertJob(newPos, job);
						double obj = neighbor.getObjective();
						double quickEval = neighbor.getMaxTardiness();
						String moveKey = "INSERT:" + job.id + ":" + i + "->" + newPos;
						Candidate cand = new Candidate(neighbor, moveKey, obj, quickEval, true);
						allCandidates.add(cand);
					}
				}
			}

			// --- Dynamic Aspiration Threshold ---
			// Base threshold is 5%. If non-improvement is high, relax the threshold (e.g., from 0.05 down to 0.01).
			double baseAspiration = 0.05;
			double dynamicAspiration = baseAspiration - ((double) nonImprovementCount / noImprovementThreshold) * 0.04;
			dynamicAspiration = Math.max(dynamicAspiration, 0.01);

			// --- Filter Candidates based on Tabu and Aspiration, and apply LTM penalty ---
			List<Candidate> validCandidates = new ArrayList<>();
			for (Candidate cand : allCandidates) {
				boolean isTabu = tabuList.containsKey(cand.moveKey) && tabuList.get(cand.moveKey) > iteration;
				// Incorporate long-term memory: add a penalty proportional to the move frequency.
				int freq = frequencyMap.getOrDefault(cand.moveKey, 0);
				double adjustedObj = cand.objectiveValue + PENALTY_COEFFICIENT * freq;
				cand.adjustedObjective = adjustedObj;

				// Allow the candidate if not tabu or if it meets the dynamic aspiration criterion.
				if (!isTabu || adjustedObj <= bestObjective * (1 - dynamicAspiration)) {
					validCandidates.add(cand);
				}
			}
			List<Candidate> candidateList = validCandidates.isEmpty() ? allCandidates : validCandidates;

			// --- Two-Phase Candidate Selection ---
			// Phase 1: sort candidates by quick evaluation (maximum tardiness).
			candidateList.sort(Comparator.comparingDouble(c -> c.quickEvaluation));
			int subsetSize = Math.max(1, candidateList.size() / 2);
			List<Candidate> candidateSubset = candidateList.subList(0, subsetSize);
			// Phase 2: select candidate with the best (lowest) adjusted objective.
			Candidate bestCandidate = Collections.min(candidateSubset, Comparator.comparingDouble(c -> c.adjustedObjective));

			// Update current solution.
			currentSolution = bestCandidate.solution;
			double currentObjective = bestCandidate.objectiveValue;

			// --- Update Long Term Memory Frequency ---
			frequencyMap.put(bestCandidate.moveKey, frequencyMap.getOrDefault(bestCandidate.moveKey, 0) + 1);

			// --- Adaptive Parameter Updates ---
			if (currentObjective < bestObjective) {
				bestSolution = new Solution(currentSolution);
				bestObjective = currentObjective;
				bestIteration = iteration;
				nonImprovementCount = 0;
				currentTabuTenure = minTabuTenure;
				currentJumpSize = minJumpSize;  // reset jump size upon improvement
			} else {
				nonImprovementCount++;
				if (nonImprovementCount >= noImprovementThreshold) {
					// Increase both tabu tenure and jump size to diversify search.
					currentTabuTenure = Math.min(currentTabuTenure + 1, maxTabuTenure);
					currentJumpSize = Math.min(currentJumpSize + 1, maxJumpSize);
					nonImprovementCount = 0;
				}
			}

			// --- Update Tabu List ---
			tabuList.put(bestCandidate.moveKey, iteration + currentTabuTenure);
			int finalIteration = iteration;
			tabuList.entrySet().removeIf(entry -> entry.getValue() <= finalIteration);
		}

		System.out.println("Best solution found at iteration: " + bestIteration);
		return bestSolution;
	}

	// Helper factorial function (to limit maxIterations).
	public static int fattoriale(int n, int maxValue) {
		if (n <= 1) return 1;
		int result = n;
		for (int i = n - 1; i >= 2; i--) {
			if (result > maxValue / i) return result;
			result *= i;
		}
		return result;
	}
}

// -------------------- Supporting Classes --------------------

// Job class remains as before.
class Job {
	int id;
	int processingTime;
	int dueDate;
	int weight;

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

// The Solution class encapsulates the schedule and its operations.
class Solution {
	List<Job> schedule;

	public Solution(List<Job> jobs) {
		this.schedule = new ArrayList<>(jobs);
	}

	// Copy constructor.
	public Solution(Solution other) {
		this.schedule = new ArrayList<>(other.schedule);
	}

	public void sortByDueDate() {
		schedule.sort(Comparator.comparingInt(j -> j.dueDate));
	}

	// Compute total weighted tardiness.
	public double getObjective() {
		double total = 0;
		int currentTime = 0;
		for (Job job : schedule) {
			currentTime += job.processingTime;
			total += job.weight * Math.max(0, currentTime - job.dueDate);
		}
		return total;
	}

	// Return maximum tardiness in the schedule.
	public double getMaxTardiness() {
		double maxTardiness = 0;
		int currentTime = 0;
		for (Job job : schedule) {
			currentTime += job.processingTime;
			double tardiness = Math.max(0, currentTime - job.dueDate);
			if (tardiness > maxTardiness) maxTardiness = tardiness;
		}
		return maxTardiness;
	}

	// Compute tardiness for the job at a given index.
	public double getJobTardiness(int index) {
		int currentTime = 0;
		for (int i = 0; i <= index; i++) {
			currentTime += schedule.get(i).processingTime;
		}
		Job job = schedule.get(index);
		return Math.max(0, currentTime - job.dueDate);
	}

	public int size() {
		return schedule.size();
	}

	public void swap(int i, int j) {
		Collections.swap(schedule, i, j);
	}

	// Remove and return job at a specified index.
	public Job removeJob(int index) {
		return schedule.remove(index);
	}

	// Insert a job at the specified index.
	public void insertJob(int index, Job job) {
		schedule.add(index, job);
	}

	// Get the job ID at a given position.
	public int getJobIdAt(int index) {
		return schedule.get(index).id;
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

// Candidate class holds a candidate solution along with its move key,
// full objective value, quick evaluation (maximum tardiness), and whether it is an insertion move.
class Candidate {
	Solution solution;
	String moveKey;   // e.g., "SWAP:1:2" or "INSERT:3:5->2"
	double objectiveValue;    // full objective value (total weighted tardiness)
	double quickEvaluation;   // quick evaluation (max tardiness)
	double adjustedObjective; // objective value plus LTM penalty
	boolean isInsertion;

	public Candidate(Solution solution, String moveKey, double objectiveValue, double quickEvaluation, boolean isInsertion) {
		this.solution = solution;
		this.moveKey = moveKey;
		this.objectiveValue = objectiveValue;
		this.quickEvaluation = quickEvaluation;
		this.isInsertion = isInsertion;
	}
}
