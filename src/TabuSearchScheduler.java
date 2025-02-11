import java.util.*;

public class TabuSearchScheduler {

	// Parameters
	public static final int MAX_FACTORIAL = 10_000;

	public static final double MIN_TENURE_FRAC = 0.005;
	public static final double INITIAL_TENURE_FRAC = 0.15;
	public static final double MAX_TENURE_FRAC = 0.75;

	public static final int NO_IMPROVEMENT_THRESHOLD_FRAC = 100;
	public static final int MIN_NO_IMPROVEMENT_THRESHOLD = 50;

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

		// Setup of parameters
		int size = jobs.size();
		int maxIterations = fattoriale(size, MAX_FACTORIAL); // maxIterations = min(n!, MAX_FCT_VALUE)

		int maxTabuTenure = (int) (jobs.size() * MAX_TENURE_FRAC); // Maximum tabu tenure
		int minTabuTenure = (int) Math.max(size * MIN_TENURE_FRAC, 1); // Minimum tabu tenure
		// Initial tabu tenure: between minTabuTenure and maxTabuTenure, prefer the lower bound to allow diversification
		int initialTabuTenure = (int) Math.max(minTabuTenure, Math.min(size * INITIAL_TENURE_FRAC, maxTabuTenure));

		int noImprovementThreshold = Math.max(maxIterations / NO_IMPROVEMENT_THRESHOLD_FRAC, MIN_NO_IMPROVEMENT_THRESHOLD);

		// Run the enhanced Tabu Search
		List<Job> bestSequence = tabuSearch(jobs, maxIterations, initialTabuTenure, minTabuTenure, maxTabuTenure, noImprovementThreshold);
	}

	public static List<Job> tabuSearch(List<Job> jobs, int maxIterations,
									   int initialTabuTenure, int minTabuTenure, int maxTabuTenure, int noImprovementThreshold) {
		// INITIAL SOLUTION (using EDD rule): sort jobs by due date
		List<Job> currentSolution = new ArrayList<>(jobs);
		currentSolution.sort(Comparator.comparingInt(job -> job.dueDate));

		List<Job> bestSolution = new ArrayList<>(currentSolution);
		double bestObjective = calculateTotalWeightedTardiness(bestSolution);

		System.out.println("EDD Initial Solution:");
		for (Job job : currentSolution) {
			System.out.print("Job" + job.id + " ");
		}
		System.out.println("\nTotal Weighted Tardiness: " + bestObjective);

		// Keep track of the best solution and its iteration
		int bestIteration = 0;

		// Tabu list: (key=move, value=iteration until which the move remains tabu)
		Map<Move, Integer> tabuList = new HashMap<>();
		int iteration = 0;

		// Adaptive tenure variables
		int currentTabuTenure = initialTabuTenure;
		int nonImprovementCount = 0;

		while (iteration < maxIterations) {
			iteration++;
			List<Candidate> candidateList = new ArrayList<>();

			// NEIGHBORHOOD GENERATION: using non-adjacent swaps
			// Also compute both the full objective and a quick evaluation (max tardiness)
			for (int i = 0; i < currentSolution.size() - 1; i++) {
				for (int j = i + 1; j < currentSolution.size(); j++) {
					List<Job> neighbor = new ArrayList<>(currentSolution);
					Collections.swap(neighbor, i, j);

					Move move = new Move(currentSolution.get(i).id, currentSolution.get(j).id);
					double neighborObjective = calculateTotalWeightedTardiness(neighbor);
					double quickEval = calculateMaxTardiness(neighbor);

					boolean isTabu = tabuList.containsKey(move) && tabuList.get(move) > iteration;

					// 3. THRESHOLD ASPIRATION:
					// Allow a move that is tabu if it improves at least 5% over the best solution.
					if (!isTabu || neighborObjective <= bestObjective * 0.95) {
						candidateList.add(new Candidate(neighbor, move, neighborObjective, quickEval));
					}
				}
			}

			// 4. CANDIDATE LIST FALLBACK:
			// If no candidate meets the criteria, ignore tabu restrictions and generate the candidate list.
			if (candidateList.isEmpty()) {
				for (int i = 0; i < currentSolution.size() - 1; i++) {
					for (int j = i + 1; j < currentSolution.size(); j++) {
						List<Job> neighbor = new ArrayList<>(currentSolution);
						Collections.swap(neighbor, i, j);

						Move move = new Move(currentSolution.get(i).id, currentSolution.get(j).id);
						double neighborObjective = calculateTotalWeightedTardiness(neighbor);
						double quickEval = calculateMaxTardiness(neighbor);
						candidateList.add(new Candidate(neighbor, move, neighborObjective, quickEval));
					}
				}
			}

			// 5. TWO-PHASE SELECTION:
			// Phase 1: Quick evaluation based on maximum tardiness.
			candidateList.sort(Comparator.comparingDouble(c -> c.quickEvaluation));
			// For example, select the best 50% of candidates according to quick eval.
			int subsetSize = Math.max(1, candidateList.size() / 2);
			List<Candidate> candidateSubset = candidateList.subList(0, subsetSize);
			// Phase 2: From the subset, select the candidate with the best (lowest) full objective value.
			Candidate bestCandidate = Collections.min(candidateSubset, Comparator.comparingDouble(c -> c.objectiveValue));

			// Update current solution with the selected move
			currentSolution = bestCandidate.sequence;
			double currentObjective = bestCandidate.objectiveValue;

			// ADAPTIVE TENURE UPDATE:
			// If we improved, record the iteration, reset the non-improvement counter, and reduce tenure.
			if (currentObjective < bestObjective) {
				bestSolution = new ArrayList<>(currentSolution);
				bestObjective = currentObjective;
				bestIteration = iteration;
				nonImprovementCount = 0;
				currentTabuTenure = minTabuTenure;
			} else {
				nonImprovementCount++;
				if (nonImprovementCount >= noImprovementThreshold) {
					// Increase the tabu tenure to force diversification.
					currentTabuTenure = Math.min(currentTabuTenure + 1, maxTabuTenure);
					nonImprovementCount = 0;
				}
			}

			// Update the Tabu List: mark the applied move as tabu for the next 'currentTabuTenure' iterations.
			tabuList.put(bestCandidate.move, iteration + currentTabuTenure);
			// Remove expired moves from the Tabu List.
			int finalIteration = iteration;
			tabuList.entrySet().removeIf(entry -> entry.getValue() <= finalIteration);
		}

		// Print out the details.
		System.out.println("Best solution found at iteration: " + bestIteration);
		System.out.println("Best Job Sequence:");
		for (Job job : bestSolution) {
			System.out.print("Job" + job.id + " ");
		}
		System.out.println("\nTotal Weighted Tardiness: " + bestObjective);

		return bestSolution;
	}

	// Full objective: total weighted tardiness.
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

	// Quick evaluation: maximum tardiness in the sequence.
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

	// A helper factorial function (used only to compute maxIterations)
	public static int fattoriale(int n, int maxValue) {
		if (n <= 1) {
			return 1;
		}

		int result = n;
		for (int i = n - 1; i >= 2; i--) {
			if (result > maxValue / i) { // Prevent fictitious overflow before multiplication
				return result;
			}
			result *= i;
		}
		return result;
	}
}

//------------------ Supporting Classes ------------------

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
}

class Move {
	int jobId1;
	int jobId2;

	public Move(int jobId1, int jobId2) {
		// Order the job ids consistently.
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

// Extended Candidate class to hold both full objective and quick evaluation measure.
class Candidate {
	List<Job> sequence;
	Move move;
	double objectiveValue;    // full objective (total weighted tardiness)
	double quickEvaluation;   // quick evaluation (max tardiness)

	public Candidate(List<Job> sequence, Move move, double objectiveValue, double quickEvaluation) {
		this.sequence = sequence;
		this.move = move;
		this.objectiveValue = objectiveValue;
		this.quickEvaluation = quickEvaluation;
	}
}