import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Stack;

class HairPin {
    private static final int k = 10;
    private static final int maxLoopLength = 50;
    private static final int maxHairpinLength = 400;
    private static final int minHairpinLength = 200;
    private static final int alpha = 4;
    private static final int beta = 9;

    private static final int NEITHER = 0;
    private static final int LEFT = 1;
    private static final int UP = 2;
    private static final int UPRIGHT = 3;


    public static int hasKmerAndReverse(String in, int begin) {
        if (in.length() - begin < maxLoopLength)
            return -1;

        String kMer = in.substring(begin, begin+10);
        int kMerEnds = begin + 10;
        String reversed = new StringBuffer(kMer).reverse().toString();
        for(int i=kMerEnds + (maxLoopLength - 2*k);i<begin + maxHairpinLength - k;i++) {
            if (i+k >= in.length())
                break;
            String s = in.substring(i, i+k);
            if (reversed.equals(s))
                return i;
        }
        return 0;
    }

    public static HashMap<String, String> kMerAndReverse(String in, int begin, int reverseBegin) {
        String kMer = in.substring(begin, begin + 10);
        String reversed = new StringBuffer(kMer).reverse().toString();
        HashMap<String, String> res = new HashMap<>();
        res.put("k-mer", kMer);
        res.put("k-mer begins", begin + "");
        res.put("reversed k-mer", reversed);
        res.put("reverse begins", reverseBegin+"");
        return res;
    }

    public static boolean hasCandidateRegion(String in, HashMap<String, String> kMerResult) {
        int candidateRegionLength = maxHairpinLength;
        int kMerBegins = Integer.parseInt(kMerResult.get("k-mer begins"));
        int reverseKmerBegins = Integer.parseInt(kMerResult.get("reverse begins"));
        int reverseKmerEnds = reverseKmerBegins + k;
        int lenSM = reverseKmerBegins - (kMerBegins + k);
        int lenSLR = (candidateRegionLength - (lenSM + 2*k))/2;
        boolean res = reverseKmerEnds + lenSLR <= in.length() -1;

        return res;
    }

    public static HashMap<String, String> cutOutCandidateRegion(String in, HashMap<String, String> kMerResult) {
        HashMap<String, String> res = new HashMap<>();

        int kMerbegins = Integer.parseInt(kMerResult.get("k-mer begins"));
        int reverseKmerBegins = Integer.parseInt(kMerResult.get("reverse begins"));
        int lenSM = reverseKmerBegins - (kMerbegins + k);
        int lenSL = (maxHairpinLength - (lenSM + 2*k))/2;
        int slBegins = kMerbegins - lenSL;
        int srBegins = reverseKmerBegins + k;
        int smBegins = slBegins + lenSL + k;
        int candidateLength = 2*lenSL + 2*k + lenSM;
        String kMer = kMerResult.get("k-mer");
        String reverseKMer = kMerResult.get("reversed k-mer");
        String sl = in.substring(slBegins, slBegins + lenSL);
        String sr = in.substring(srBegins, srBegins + lenSL);
        String srReversed = new StringBuffer(sr).reverse().toString();
        String sm = in.substring(smBegins, smBegins+lenSM);
        String smReversed = new StringBuffer(sm).reverse().toString();
        res.put("k-mer", kMer);
        res.put("reverse k-mer", reverseKMer);
        res.put("SL begins", slBegins+"");
        res.put("k-mer begins", kMerbegins+"");
        res.put("SL", sl);
        res.put("SM", sm);
        res.put("SR reversed", srReversed);
        res.put("SM reversed", smReversed);
        res.put("candidate length", candidateLength+"");
        return res;
    }

    public static int[][] findLCSArrows(String a, String b) {
        int lenA = a.length();
        int lenB = b.length();
        int[][] numbers = new int[lenA+1][lenB+1];
        int[][] arrows = new int[lenA+1][lenB+1];

        for(int i=0;i<=lenA;i++) {
            numbers[i][0] = 0;
            arrows[i][0] = NEITHER;
        }

        for(int j=0;j<=lenB;j++) {
            numbers[0][j] = 0;
            arrows[0][j] = NEITHER;
        }

        for(int j=1;j<=lenB;j++) {
            for(int i=1;i<=lenA;i++) {
                if (a.charAt(i-1)==b.charAt(j-1)) {
                    numbers[i][j] = numbers[i-1][j-1] + 1;
                    arrows[i][j] = UPRIGHT;
                } else if (numbers[i-1][j] >= numbers[i][j-1]) {
                    numbers[i][j] = numbers[i-1][j];
                    arrows[i][j]=UP;
                } else {
                    numbers[i][j] = numbers[i][j-1];
                    arrows[i][j] = LEFT;
                }
            }
        }
        return arrows;
    }

    public static void printLCS(String a, int[][] arrows, int i, int j) {
        if (i==1 || j==1)
            return;
        if(arrows[i-1][j-1]==UPRIGHT) {
            printLCS(a, arrows, i-1, j-1);
            System.out.print(a.charAt(i-2));
        } else if (arrows[i-1][j-1]==UP) {
            printLCS(a, arrows, i-1, j);
        }
        else
            printLCS(a, arrows, i, j-1);
    }

    public static String[] align(int[][] arrows, String a, String b) {
        int row = a.length()+1;
        int col = b.length()+1;
        Stack<Character> stack = new Stack<>();

        // row
        while(row > 0 && col > 0) { // || ? &&
            if (arrows[row-1][col-1] == UPRIGHT) {
                stack.push(a.charAt(row-2));
                row-=1;
                col-=1;
            } else if (arrows[row-1][col-1] == UP) {
                stack.push(a.charAt(row-2));
                row-=1;
            } else if (arrows[row-1][col-1] == LEFT) {
                stack.push('-');
                col-=1;
            } else if (arrows[row-1][col-1] == NEITHER) {
                // one of cases (0, 0) , (1, 0), (0, 1)
                if (col == 1 && row == 1)
                    break;
                else if (col == 1 && row > 1) {
                    stack.push(a.charAt(row - 2));
                    row -= 1;
                } else if ( row == 1 && col > 1) {
                    stack.push('-');
                    col -= 1;
                }

            }
        }

        StringBuilder sbl = new StringBuilder();
        while(!stack.empty()) {
            sbl.append(stack.pop());
        }
        String sl = sbl.toString();
        stack.removeAllElements();

        row = a.length() + 1;
        col = b.length() + 1;
        // col
        while(row > 0 && col > 0) {
            if (arrows[row-1][col-1] == UPRIGHT) {
                stack.push(b.charAt(col-2));
                row-=1;
                col-=1;
            } else if (arrows[row-1][col-1]==UP) {
                stack.push('-');
                row-=1;
            } else if (arrows[row-1][col-1]==LEFT) {
                stack.push(b.charAt(col-2));
                col-=1;
            } else if (arrows[row-1][col-1] == NEITHER) {
                if (col == 1 && row == 1)
                    break;
                else if (col == 1 && row > 1) {
                    stack.push('-');
                    row -= 1;
                } else if (row == 1 && col > 1) {
                    stack.push(b.charAt(col-2));
                    col -= 1;
                }
            }
        }

        StringBuilder sbr = new StringBuilder();
        while(!stack.empty()) {
            sbr.append(stack.pop());
        }
        String sr = sbr.toString();
        String[] res = new String[2];
        res[0] = sl;
        res[1] = sr;
        return res;
    }

    public static HashMap<String, String> folding(HashMap<String, String> cutOutted) {
        String kMer = cutOutted.get("k-mer");
        String sl = cutOutted.get("SL");
        String srReversed = cutOutted.get("SR reversed");
        String sm = cutOutted.get("SM");
        String smReversed = cutOutted.get("SM reversed");
        int[][] sLRreverseArrows = findLCSArrows(sl, srReversed);
        int[][] sMandsMreversedArrows = findLCSArrows(sm, smReversed);
        String[] sLRaligned = align(sLRreverseArrows, sl, srReversed);
        String[] sMaligned = align(sMandsMreversedArrows, sm, smReversed);

        String upper = sLRaligned[0] + kMer + sMaligned[0];
        String lower = sLRaligned[1] + kMer + sMaligned[1];
        HashMap<String, String> res = new HashMap<>();
        res.put("upper", upper);
        res.put("lower", lower);
        res.put("k-mer on lower begins", sLRaligned[0].length()+"");
        res.put("k-mer", kMer);
        res.put("k-mer reversed", new StringBuilder(kMer).reverse().toString());
        return res;
    }

    public static HashMap<String, String> candidateHairpin(HashMap<String, String> folded) {
        String upper = folded.get("upper");
        String lower = folded.get("lower");
        int kMerOnUpperBegins = Integer.parseInt(folded.get("k-mer on lower begins"));
        int kMerOnUpperEnds = kMerOnUpperBegins + k;

        int leftStop = kMerOnUpperBegins;
        int leftPrev = kMerOnUpperBegins;
        int leftLastMatch = kMerOnUpperBegins;

        boolean matchPosition = upper.charAt(leftStop) != '-' && lower.charAt(leftStop) != '-';
        boolean prevIndelPosition = upper.charAt(leftPrev) == '-' || lower.charAt(leftPrev) == '-';
        boolean indelNumberCondition = indelNumberCondition(upper, lower, kMerOnUpperBegins);

        while(!matchPosition || !prevIndelPosition || !indelNumberCondition) {
            if (matchPosition)
                leftLastMatch = leftStop;
            leftPrev = leftStop;
            leftStop -= 1;

            if (leftStop < 0) {
                break;
            }

            matchPosition = upper.charAt(leftStop) != '-' && lower.charAt(leftStop) != '-';
            prevIndelPosition = upper.charAt(leftPrev) == '-' || lower.charAt(leftPrev) == '-';
            indelNumberCondition = indelNumberCondition(upper, lower, leftStop);
        }

        int rightStop = kMerOnUpperEnds;
        int rightPrev = kMerOnUpperEnds;
        int rightLastMatch = kMerOnUpperEnds;
        boolean rMatchPosition = upper.charAt(rightStop) != '-' && lower.charAt(rightStop) != '-';
        boolean rPrevIndelPosition = upper.charAt(rightPrev) == '-' || lower.charAt(rightPrev) == '-';
        boolean rIndelNumberCondition = indelNumberCondition(upper, lower, kMerOnUpperEnds);

        while(!rMatchPosition || !rPrevIndelPosition || !rIndelNumberCondition) {
            if (rMatchPosition) {
                rightLastMatch = rightStop;
            }
            rightPrev = rightStop;
            rightStop += 1;


            if (rightStop > upper.length() - 1) {
                break;
            }

            rMatchPosition = upper.charAt(rightStop) != '-' && lower.charAt(rightStop) != '-';
            rPrevIndelPosition = upper.charAt(rightPrev) == '-' || lower.charAt(rightPrev) == '-';
            rIndelNumberCondition = indelNumberCondition(upper, lower, rightStop);
        }

        boolean stopIdxInArms = (leftStop == kMerOnUpperBegins) || (rightStop == kMerOnUpperEnds);

        String leftArm = upper.substring(leftLastMatch, rightLastMatch + 1);
        leftArm = leftArm.replaceAll("-", "");
        String rightArm = lower.substring(leftLastMatch, rightLastMatch + 1);
        rightArm = rightArm.replaceAll("-", "");
        String rightReversed = new StringBuilder(rightArm).reverse().toString();
        String loop = upper.substring(rightLastMatch + 1, upper.length() - (rightLastMatch - kMerOnUpperEnds)-1);
        loop = loop.replaceAll("-", "");
        String stopIdxArmsCondition;
        if (stopIdxInArms) {
            stopIdxArmsCondition = "invalid";
        } else {
            stopIdxArmsCondition = "valid";
        }
        HashMap<String, String> res = new HashMap<String, String>();
        res.put("left", leftArm);
        res.put("right", rightArm);
        res.put("loop", loop);
        res.put("hairpin", leftArm + loop + rightReversed);
        res.put("stopIdxArmsCondition", stopIdxArmsCondition);

        return res;

    }

    private static boolean indelNumberCondition(String upper, String lower, int current) {

        int indelCount=0;

        for(int i=current-1; i>=current-beta/2;i--) {
            if (i >= 0) {
                if (upper.charAt(i)=='-' || lower.charAt(i)=='-')
                    indelCount += 1;
            }
        }

        for (int i=current+1; i<=current+beta/2;i++) {
            if (i > 0) {
                if (upper.charAt(i) == '-' || lower.charAt(i)=='-')
                    indelCount += 1;
            }
        }

        return indelCount >= alpha;
    }

    public static boolean verifyHairpin(String leftArm, String rightArm, String loop) {
        String hairpin = leftArm + loop + rightArm;
        boolean loopLengthCondition = loop.length() >= 0 && loop.length() <= maxLoopLength;
        boolean hairpinLengthCondition = hairpin.length() >= minHairpinLength && hairpin.length() <= maxHairpinLength;

        return loopLengthCondition && hairpinLengthCondition;

    }

}

public class HW3 {
    public static void main(String[] args) throws IOException {
        FileInputStream in = new FileInputStream(args[0]);
        Scanner sc = new Scanner(in);
        String title = sc.nextLine();
        String sequence = sc.nextLine();

        int begin = 0;
        while(true) {
            int hasKmerAndReverse = HairPin.hasKmerAndReverse(sequence, begin);
            if (hasKmerAndReverse == 0) {
                begin += 1;
                continue;
            } else if (hasKmerAndReverse == -1 ) {
                break;
            } else {
                HashMap<String, String> kMerFound = HairPin.kMerAndReverse(sequence, begin, hasKmerAndReverse);
                HashMap<String, String> cutOut;
                if (HairPin.hasCandidateRegion(sequence, kMerFound)) {
                    cutOut = HairPin.cutOutCandidateRegion(sequence, kMerFound);
                } else {
                    break;
                }

                HashMap<String, String> hairpinCandidate = HairPin.candidateHairpin(HairPin.folding(cutOut));
                boolean verified = HairPin.verifyHairpin(hairpinCandidate.get("left"), hairpinCandidate.get("right"), hairpinCandidate.get("loop"));
                if (verified) {
                    System.out.println(hairpinCandidate.get("hairpin"));
                    HairPin.printLCS(hairpinCandidate.get("left"),
                            HairPin.findLCSArrows(hairpinCandidate.get("left"), hairpinCandidate.get("right")),
                            hairpinCandidate.get("left").length() + 1  ,
                            hairpinCandidate.get("right").length() + 1);
                    System.out.println();
                    System.out.println(hairpinCandidate.get("loop"));
                    System.out.println();
                }
                begin = Integer.parseInt(cutOut.get("SL begins")) + Integer.parseInt(cutOut.get("candidate length"));
            }
        }
    }
}

