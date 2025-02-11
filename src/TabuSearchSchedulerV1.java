import java.util.*;

public class TabuSearchSchedulerV1 {

    public static long fattoriale(int n) {
        if (n <= 1) {
            return 1;
        } else {
            for (int i = n - 1; i > 1; i--) {
                n *= i;
                if (n > 10_000) {
                    return n;
                }
            }
        }
        return n;
    }

    public static void main(String[] args) {
        int w=1;
        // Define jobs with processing time (p), due date (d), and weight (w)
        List<Job> jobs = Arrays.asList(
                new Job(1, 6, 9, w),
                new Job(2, 4, 12, w),
                new Job(3, 8, 15, w),
                new Job(4, 2, 8, w),
                new Job(5, 10, 20, w),
                new Job(6, 3, 22, w)
        );

        // min(sqrt(n!), 10_000)
        int maxIterations = (int) Math.min(Math.sqrt(fattoriale(jobs.size())), 10_000); // Maximum number of iterations
        int tabuFactor = 20;      // Tabu tenure factor
        int tabuTenure = Math.min(maxIterations/tabuFactor, jobs.size()-1);      // Tabu tenure
        System.out.println("Tabu tenure: " + tabuTenure);

        // Run Tabu Search
        List<Job> bestSequence = tabuSearch(jobs, maxIterations, tabuTenure);

        // Print the best sequence found
        System.out.println("Best Job Sequence:");
        for (Job job : bestSequence) {
            System.out.print("Job" + job.id + " ");
        }

        // Calculate and print total weighted tardiness
        double totalWeightedTardiness = calculateTotalWeightedTardiness(bestSequence);
        System.out.println("\nTotal Weighted Tardiness: " + totalWeightedTardiness);
    }

    public static List<Job> tabuSearch(List<Job> jobs, int maxIterations, int tabuTenure) {
        // Initial solution using EDD rule (Earliest Due Date)
        List<Job> currentSolution = new ArrayList<>(jobs);
        currentSolution.sort(Comparator.comparingInt(job -> job.dueDate));

        List<Job> bestSolution = new ArrayList<>(currentSolution);
        double bestObjective = calculateTotalWeightedTardiness(bestSolution);

        Map<Move, Integer> tabuList = new HashMap<>();
        int iteration = 0;

        while (iteration < maxIterations) {
            iteration++;

            List<Candidate> candidateList = new ArrayList<>();

            // Generate neighborhood by swapping adjacent jobs
            for (int i = 0; i < currentSolution.size() - 1; i++) {
                List<Job> neighbor = new ArrayList<>(currentSolution);
                Collections.swap(neighbor, i, i + 1);

                Move move = new Move(currentSolution.get(i).id, currentSolution.get(i + 1).id);
                double neighborObjective = calculateTotalWeightedTardiness(neighbor);

                boolean isTabu = tabuList.containsKey(move) && tabuList.get(move) > iteration;

                // Aspiration criterion
                if (!isTabu || neighborObjective < bestObjective) {
                    candidateList.add(new Candidate(neighbor, move, neighborObjective, 0.0));
                }
            }

            if (candidateList.isEmpty()) {
                break; // No feasible moves
            }

            // Select the best candidate
            Candidate bestCandidate = Collections.min(candidateList, Comparator.comparingDouble(c -> c.objectiveValue));

            currentSolution = bestCandidate.sequence;
            double currentObjective = bestCandidate.objectiveValue;

            // Update best solution
            if (currentObjective < bestObjective) {
                bestSolution = new ArrayList<>(currentSolution);
                bestObjective = currentObjective;
            }

            // Update Tabu List
            tabuList.put(bestCandidate.move, iteration + tabuTenure);

            // Remove expired moves from Tabu List
            int finalIteration = iteration;
            tabuList.entrySet().removeIf(entry -> entry.getValue() <= finalIteration);
        }

        return bestSolution;
    }

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
}
