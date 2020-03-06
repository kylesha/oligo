
/**
 * @author Ky Sha
 *
 A class dedicated to manipulating/searching Oligo objects. All coordinates are ZERO-based
 *
 * METHODS ON REQUEST replace() mutateINS() - insertional mutation mutateDEL() - deletion mutation mutateTranspose() - transpositional insertion label5Prime() - adds
 * adduct (i.e. biotin) to 5 prime end label3Prime() - adds adduct to 3 primer end explicitMode() - show 5'/3' labels of specified oligo constructors strip away
 * punctuation marks (including space) embedded in input oligo -------------------------------------------------------------------------------------------------
 */
import com.google.common.base.Preconditions;
import java.util.*;
import java.io.*;
import com.google.common.primitives.*;
import java.util.regex.*;

public class Oligo implements Serializable
{
	private final String oligo;
	private final String mnf = "MatchNotFound"; //exception
	private final String soob = "StartIndexOutOfBounds"; //exception
	private final String eoob = "EndIndexOutOfBounds"; //exception
	private char ignoredChar = 'n'; //default character to ignore in oligo
	private final int oligo_length;
	private static final String OLIGO_ALPHABET_REGEX = "[acgtnACGTN]+";

//============================| CONSTRUCTORS |================================//

	/**
	 * Constructor: Constructs an oligo object of arbitrary length
	 *
	 * @param oligo (type: String, StringBuilder, StringBuffer, void)
	 */
	public Oligo(String oligo)
	{
		Pattern pat = Pattern.compile(OLIGO_ALPHABET_REGEX);
		Matcher mat = pat.matcher(oligo);

		Preconditions.checkArgument(mat.matches(), "Oligo object cannot be created. Input Oligo object [%s] contains invalid characters.", oligo);
		this.oligo = oligo;
		oligo_length = oligo.length();
	}


	public Oligo(StringBuilder oligo)
	{
		this.oligo = oligo.toString();
		oligo_length = oligo.length();
	}


	public Oligo(StringBuffer oligo)
	{
		this.oligo = oligo.toString();
		oligo_length = oligo.length();
	}


	public Oligo(Oligo oligo)
	{
		this.oligo = oligo.toString();
		oligo_length = oligo.length();
	}


	public Oligo()
	{
		oligo = "";
		oligo_length = oligo.length();
	}


//=========================| OVER RIDDEN METHODS |============================//
	/**
	 * Converts oligo to character array
	 *
	 * @return CharArray representation of oligo
	 */
	public char[] toCharArray()
	{
		return oligo.toCharArray();
	}


	/**
	 * Converts all characters in oligo to lowercase characters
	 *
	 * @return Oligo object - oligo in lowercase
	 */
	public Oligo toLowerCase()
	{
		return new Oligo(oligo.toLowerCase());
	}


	/**
	 * Converts oligo to String object
	 *
	 * @return String object representation of oligo
	 */
	@Override
	public String toString()
	{
		return oligo;
	}


	/**
	 * Converts oligo to StringBuilder object
	 *
	 * @return StringBuilder object representation of oligo
	 */
	public StringBuilder toStringBuilder()
	{
		return new StringBuilder(oligo);
	}


	public Oligo toUpperCase()
	{
		return new Oligo(oligo.toUpperCase());
	}


//=================================| METHODS |================================//
	/**
	 * Returns the reverse complement of the input oligo
	 *
	 * @return
	 */
	public Oligo antiparallel()
	{
		Oligo tempOligo = new Oligo(oligo);
		Oligo complement = new Oligo(tempOligo.complement());
		return complement.reverse();
	}


	public TreeMap<String, Integer> baseContent()
	{
		int gCount = 0;
		int aCount = 0;
		int tCount = 0;
		int cCount = 0;
		char[] source = oligo.toCharArray();

		for(int i = 0; i <= oligo_length - 1; i++)
		{
			if(source[i] == 'G' || source[i] == 'g')
				gCount++;
			else if(source[i] == 'A' || source[i] == 'a')
				aCount++;
			else if(source[i] == 'T' || source[i] == 't')
				tCount++;
			else if(source[i] == 'C' || source[i] == 'c')
				cCount++;
		}
		TreeMap<String, Integer> output = new TreeMap<String, Integer>();

		output.put("A", aCount);
		output.put("C", cCount);
		output.put("G", gCount);
		output.put("T", tCount);

		return output;
	}


	/**
	 * Returns the complement of the input oligo object
	 *
	 * @return Oligo object
	 *
	 */
	public Oligo complement()
	{
		String tempOligo = "";
		char[] source = oligo.toCharArray();

		for(int i = 0; i <= oligo_length - 1; i++)
		{
			if(source[i] == 'a')
				tempOligo += 't';
			else if(source[i] == 'A')
				tempOligo += 'T';
			else if(source[i] == 'c')
				tempOligo += 'g';
			else if(source[i] == 'C')
				tempOligo += 'G';
			else if(source[i] == 'g')
				tempOligo += 'c';
			else if(source[i] == 'G')
				tempOligo += 'C';
			else if(source[i] == 't')
				tempOligo += 'a';
			else if(source[i] == 'T')
				tempOligo += 'A';
			else
				tempOligo += source[i];
		}
		return new Oligo(tempOligo);
	}


	public boolean contains(String s)
	{
		return oligo.contains(s);
	}


	public boolean contains(Oligo query)
	{
		return oligo.contains(query.toString());
	}


	/**
	 * Counts the specified character in oligo object
	 *
	 * @param inputChar character to be counter
	 * @return int represent number of character counted
	 */
	public int countChar(char inputChar)
	{
		char ch = inputChar;
		int count = 0;
		for(int i = 0; i <= oligo_length - 1; i++)
		{
			if(oligo.charAt(i) == ch || oligo.charAt(i) == ch)
				count++;
		}
		return count;
	}


	/**
	 * Progressively cuts oligo at indicated coordinates and returns list of resulting cut fragments in order of increasing coordinate Largest input cut site cannot be
	 * larger than coordinate of last bp. For example, if oligo has 100 base pairs, max bp cannot be greater than 99
	 *
	 * @param bp
	 * @return
	 * @throws OligoException
	 */
	public List<Oligo> cutAt(int... bp) throws OligoException
	{
		Arrays.sort(bp);
		List<Oligo> cutFrags = new LinkedList<Oligo>();

		int max = Ints.max(bp);
		if(max <= this.oligo_length - 1)
		{
			List<Interval> cutCoordinates = new LinkedList<Interval>();

			cutCoordinates.add(new Interval(0, bp[0])); //first cut coordinate

			for(int i = 0; i + 1 < bp.length; i++)
				cutCoordinates.add(new Interval(bp[i] + 1, bp[i + 1]));

			cutCoordinates.add(new Interval(max + 1, this.oligo_length - 1)); //last cut coordinate

			for(Interval i : cutCoordinates)
				cutFrags.add(extractSequence(i.start(), i.stop()));
		}
		else
		{
			System.out.println("ERROR: [" + Thread.currentThread().getStackTrace()[1].getClassName() + "." + Thread.currentThread().getStackTrace()[1].getMethodName() + "]: Cut coordinate larger than length of oligo to be cut. Thread aborted.");
			System.exit(0);
		}
		return cutFrags;
	}


	public static String decodeLongToDNA(Long input)
	{
		char[] cha = input.toString().toCharArray();
		String dna = "";

		for(char ch : cha)
		{
			if(ch == '1')
				dna += 'a';
			else if(ch == '2')
				dna += 'c';
			else if(ch == '3')
				dna += 'g';
			else if(ch == '4')
				dna += 't';
		}
		return dna;
	}


	public static long encodeToLong(String str)
	{
		char[] cha = str.toLowerCase().toCharArray();
		String strInt = "";

		for(char ch : cha)
		{
			if(ch == 'a')
				strInt += '1';
			else if(ch == 'c')
				strInt += '2';
			else if(ch == 'g')
				strInt += '3';
			else if(ch == 't')
				strInt += '4';
		}

		long dna = Long.parseLong(strInt);
		if(dna > Long.MAX_VALUE)
		{
			System.out.println("ERROR: [" + Thread.currentThread().getStackTrace()[1].getClassName() + "." + Thread.currentThread().getStackTrace()[1].getMethodName() + "]: Long DNA value exceeds Long.MAXIMUM_VALUE. Thread aborted.");
			System.exit(0);
		}
		return dna;
	}


	/**
	 * Excises everthing left of query, including the query sequence
	 *
	 * @param query search key
	 * @param mismatches max number of allowed mismatches
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseLeftFrom(Oligo query, int mismatches) throws OligoException
	{
		if(isFuzzyMatch(query, mismatches))
		{
			int index = getFirstMatchCoordinate(query, mismatches);
			return new Oligo(oligo.substring(index + query.length()));
		}
		else
			throw new OligoException(mnf, "exciseLeftFrom()"); // exhausted all possibilities, no matches found
	} // end exciseLeftFrom()


	/**
	 * Excises everthing left of query, including the query sequence
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseLeftFrom(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
	{
		if(isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength))
		{
			Oligo swKey = getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
			return new Oligo(exciseLeftFrom(swKey, mismatches));
		}
		else
			throw new OligoException(mnf, "exciseLeftFrom()"); // exhausted all possibilities, no matches found
	} // end exciseLeftFrom()


	/**
	 * Excises everthing left of query
	 *
	 * @param query search key
	 * @param mismatches max number of allowed mismatches
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseLeftOf(Oligo query, int mismatches) throws OligoException
	{
		if(isFuzzyMatch(query, mismatches))
		{
			int index = getFirstMatchCoordinate(query, mismatches);
			return new Oligo(oligo.substring(index));
		}
		else
			throw new OligoException(mnf, "exciseLeftOf()"); // exhausted all possibilities, no matches found
	}//end exciseLeftOf()


	/**
	 * Excises everthing left of query
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseLeftOf(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
	{
		int index = getFirstMatchCoordinate(query, mismatches, ins, del, minKeyLength);
		if(index != -1)
			return new Oligo(oligo.substring(index));
		else
			throw new OligoException(mnf, "exciseLeftOf()");
	}//end exciseLeftOf()


	/**
	 * Excises everthing right of query, including the query sequence
	 *
	 * @param query search key
	 * @param mismatches max number of allowed mismatches
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseRightFrom(Oligo query, int mismatches) throws OligoException
	{
		if(isFuzzyMatch(query, mismatches))
		{
			int index = getFirstMatchCoordinate(query, mismatches);
			return new Oligo(oligo.substring(0, index));
		}
		else
			throw new OligoException(mnf, "exciseRightFrom()"); // exhausted all possibilities, no matches found
	} //end xciseRightFrom()


	/**
	 * Excises everthing right of query, including the query sequence
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseRightFrom(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
	{
		int index = getFirstMatchCoordinate(query, mismatches, ins, del, minKeyLength);
		if(index != -1)
			return new Oligo(oligo.substring(0, index));
		else
			throw new OligoException(mnf, "exciseRightFrom()");
	}//end xciseRightFrom()


	/**
	 * Excises everthing right of query
	 *
	 * @param query search key
	 * @param mismatches max number of allowed mismatches
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseRightOf(Oligo query, int mismatches) throws OligoException
	{
		if(isFuzzyMatch(query, mismatches))
		{
			int index = getFirstMatchCoordinate(query, mismatches);
			return new Oligo(oligo.substring(0, index + query.length()));
		}
		else
			throw new OligoException(mnf, "exciseRightOf()"); // exhausted all possibilities, no matches found
	} //end exciseRightOf()


	/**
	 * Excises everthing right of query
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo exciseRightOf(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
	{
		if(isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength))
		{
			int index = getFirstMatchCoordinate(query, mismatches, ins, del, minKeyLength);
			Oligo swKey = getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
			return new Oligo(oligo.substring(0, index + swKey.length()));
		}
		else
			throw new OligoException(mnf, "exciseRightOf()");
	}//end exciseRightOf()


	/**
	 * Extracts the first best match sequence from oligo, based on fuzzy-match algorithm
	 *
	 * @param query the sequence to extract
	 * @param mismatches the maximum number of allowed mismatches
	 * @return [Oligo object] Resulting oligo sequence after excision; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo extractSequence(Oligo query, int mismatches) throws OligoException
	{
		if(isFuzzyMatch(query, mismatches))
		{
			if(getFirstMatchCoordinate(query, mismatches) != -1)
			{
				int start = getFirstMatchCoordinate(query, mismatches);
				int stop = start + query.length();
				return new Oligo(oligo.substring(start, stop));
			}
			else
				throw new OligoException(mnf, "extractSequence()");
		}
		else
			throw new OligoException(mnf, "extractSequence()");
	} //end extractSequence()


	public Oligo extractSequence(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
	{
		Oligo swKey = new Oligo();
		try
		{
			swKey = getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
		}
		catch(OligoException error)
		{
			error.getStackTrace();
		}
		return extractSequence(swKey, mismatches);
	} //end extractSequence()


	/**
	 * Extracts a sequence based on given start and end indices, inclusively
	 *
	 * @param start start index
	 * @param end end index
	 * @return [Oligo object] Extracted source sequence; otherwise throws OligoException
	 * @throws OligoException
	 */
	public Oligo extractSequence(int start, int end) throws OligoException
	{
		String s = new String();
		boolean START_INDEX_OUT_OF_BOUND = start < 0;
		boolean END_INDEX_OUT_OF_BOUND = end > oligo.length() - 1;

		if(start >= 0 && end <= oligo.length() - 1)
			s = oligo.substring(start, end + 1);
		else if(START_INDEX_OUT_OF_BOUND ^ END_INDEX_OUT_OF_BOUND)
		{
			if(START_INDEX_OUT_OF_BOUND)
				throw new OligoException(soob, "extractSequence()");
			else
				throw new OligoException(eoob, "extractSequence()");
		}
		else if(START_INDEX_OUT_OF_BOUND && END_INDEX_OUT_OF_BOUND)
			throw new OligoException("(Start&End)IndexOutofBounds", "extractSequence()");
		return new Oligo(s);
	} //end extractSequence()


	/**
	 * Generates a random base
	 *
	 * @return char
	 */
	public static char generateRandomBase()
	{
		char ch = ' ';
		int base = (int) (Math.random() * 100);

		if(base >= 0 && base <= 24)
			ch = 'a';
		else if(base >= 25 && base <= 49)
			ch = 'c';
		else if(base >= 50 && base <= 74)
			ch = 'g';
		else if(base >= 75 && base <= 99)
			ch = 't';
		return ch;
	} //end generateRandomBase()


	/**
	 * Generate a random oligo of a specified length
	 *
	 * @param length length of random oligo
	 * @return Oligo object
	 */
	public static Oligo generateRandomOligo(int length)
	{
		String randomOligo = new String();

		for(int i = 0; i <= length - 1; i++)
		{
			int base = (int) (Math.random() * 100);

			if(base >= 0 && base <= 24)
				randomOligo += 'a';
			else if(base >= 25 && base <= 49)
				randomOligo += 'c';
			else if(base >= 50 && base <= 74)
				randomOligo += 'g';
			else if(base >= 75 && base <= 99)
				randomOligo += 't';
		}
		return new Oligo(randomOligo);
	} //end generateRandomOligo()


	/**
	 * Returns an ArrayList containing the indices of all occurences of the query
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @return ArrayList object (ArrayList is empty if zero matches are found)
	 */
	@SuppressWarnings("unchecked")
	public ArrayList<Integer> getAllMatchCoordinates(Oligo query, int mismatches)
	{
		Oligo searchWindow;
		ArrayList hits = new ArrayList();
		final int QUERY_LENGTH = query.length();

		for(int i = 0; i + QUERY_LENGTH <= oligo_length; i++)
		{
			searchWindow = new Oligo(oligo.substring(i, i + QUERY_LENGTH));
			if(searchWindow.isFuzzyMatch(query, mismatches))
				hits.add(i);
		}
		return hits;
	}//end getAllMatchCoordinates()


	/**
	 * Returns an ArrayList containing the indices of all occurences of the query
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return ArrayList object (ArrayList is empty if zero matches are found)
	 * @see smithWaterman()
	 * @see isFuzzSmithWatermanMatch()
	 */
	@SuppressWarnings("unchecked")
	public ArrayList<Integer> getAllMatchCoordinates(Oligo query, int mismatches, int ins, int del, int minKeyLength)
	{
		Oligo searchWindow;
		Oligo swKey;
		ArrayList hits = new ArrayList();
		final int QUERY_LENGTH = query.length(); //net query length must include allowed number of inserts in search window

		for(int i = 0; i + QUERY_LENGTH + ins <= oligo_length; i++)
		{
			try
			{
				searchWindow = new Oligo(oligo.substring(i, i + QUERY_LENGTH + ins));   //searchWindow must compensate for number of allowed inserts
				swKey = searchWindow.smithWaterman(query, ins, del);
				searchWindow = new Oligo(oligo.substring(i, i + swKey.length()));       //narrow searchWindow to length of potential key
			}
			catch(OligoException e)
			{
				continue;
			}
			if(searchWindow.isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength) && searchWindow.isFuzzyMatch(swKey, mismatches))
				hits.add(i);
		}
		return hits;
	}//end getAllMatchCoordinates()


	/**
	 * Returns the index of the first occurence of query
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @return int index if found; otherwise returns -1
	 */
	public int getFirstMatchCoordinate(Oligo query, int mismatches)
	{
		ArrayList list = getAllMatchCoordinates(query, mismatches);
		return (list.size() > 0) ? Integer.parseInt(list.get(0).toString()) : -1;
	}//end getFirstMatchCoordinate()


	/**
	 * Returns the index of the first occurence of query
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return int index if found; otherwise returns -1
	 */
	public int getFirstMatchCoordinate(Oligo query, int mismatches, int ins, int del, int minKeyLength)
	{
		ArrayList list = getAllMatchCoordinates(query, mismatches, ins, del, minKeyLength);
		return (list.size() > 0) ? Integer.parseInt(list.get(0).toString()) : -1;
	}//end getFirstMatchCoordinate()


	/**
	 * Returns the found key (modified query) based on the given parameters. This method is based on isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int
	 * del, int minKeyLength) and simply returns the key equivalent.
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return [Oligo object] Found key
	 * @throws OligoException
	 * @see also isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int del, int minKeyLength)
	 */
	public Oligo getFuzzySWkey(Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
	{
		if(isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength))
			return smithWaterman(query, ins, del);
		else
			throw new OligoException(mnf, "getFuzzySWkey()");
	}//end getFuzzySWkey()


	/**
	 * Returns the currenly set character to ignore
	 *
	 * @return char character
	 */
	public char getIgnoredChar()
	{
		return ignoredChar;
	}


	/**
	 * Returns the start index of the last occurence of a match. If only one match found then result is same is getFirstMatchCoordinate()
	 *
	 * @param query The query to search for
	 * @param mismatches The maximum number of allowed mismatches
	 * @return int representing the found index; otherwise returns -1
	 */
	public int getLastMatchCoordinate(Oligo query, int mismatches)
	{
		ArrayList list = getAllMatchCoordinates(query, mismatches);
		return (list.size() > 0) ? Integer.parseInt(list.get(list.size() - 1).toString()) : -1;
	}


	/**
	 * Returns the start index of the last occurence of a match. If only one match found then result is same is getFirstMatchCoordinate()
	 *
	 * @param query query sequence
	 * @param mismatches maximum number of allowedMismatches
	 * @param ins maximum allowed number of inserts in SOURCE sequence
	 * @param del maximum allowed number of deletes in SOURCE sequence
	 * @param minKeyLength minimum length of found key in order for search to be considered successful
	 * @return int representing the found index; otherwise returns -1
	 */
	public int getLastMatchCoordinate(Oligo query, int mismatches, int ins, int del, int minKeyLength)
	{
		ArrayList list = getAllMatchCoordinates(query, mismatches, ins, del, minKeyLength);
		return (list.size() > 0) ? Integer.parseInt(list.get(list.size() - 1).toString()) : -1;
	}


	/**
	 * Inserts the specified insert immediately BEFORE the specified index
	 *
	 * @param insert sequence to be inserted
	 * @param index The index where insertion is to occur
	 * @return Oligo object with inserted sequence
	 */
	public Oligo insert(Oligo insert, int index)
	{
		String insertSeq = insert.toString();
		StringBuilder targetOligo = new StringBuilder(oligo);
		return new Oligo(targetOligo.insert(index, insertSeq));
	}


	/**
	 * Returns true if the input sequence contains only 'A/a', 'C/c', 'G/g', 'T/t' and the ignoredChar; Otherwise, returns false
	 *
	 * @return boolean
	 */
	public boolean isDNA()
	{
		char[] s = oligo.toCharArray();
		for(char base : s)
		{
			if(base != 'a' && base != 'c' && base != 'g' && base != 't' && base != ignoredChar)
				return false;
		}
		return true;
	}


	/**
	 * Determines whether an oligo has zero length
	 *
	 * @return boolean
	 */
	public boolean isEmpty()
	{
		return oligo_length == 0;
	}


	/**
	 * Determines whether query sequence is contained in the oligo sequence, given the allowed mismatches. Will ignore N's. Method is not case sensitive
	 *
	 * @param inputQuery The query to search for
	 * @param mismatches
	 * @return boolean
	 */
	public boolean isFuzzyMatch(Oligo inputQuery, int mismatches)
	{
		char[] query = inputQuery.toUpperCase().toCharArray();
		char[] source = oligo.toUpperCase().toCharArray();
		int ntMatches = 0; //ntMatches = number of nucleotide-nucleotide (GC or AT) matches
		int ignoredCharMatches = 0; //ignoredCharMatches = number of nucleotide-ignoredChar (i.e. 'N') matches
		final int QUERY_LENGTH = query.length;

		for(int i = 0; i + QUERY_LENGTH <= oligo_length; i++)
		{
			for(int j = 0; j <= QUERY_LENGTH - 1; j++)
			{
				boolean nnMatch = (source[i + j] == query[j]);
				boolean niMatch = (source[i + j] == ignoredChar || query[j] == ignoredChar);

				if(nnMatch && !niMatch) //count only nucleotide-nucleotide matches
					ntMatches++;
				else if(niMatch)  //count only nucleotide-ignoredChar matches
					ignoredCharMatches++;
			} //inner for() loop

			if(ntMatches >= QUERY_LENGTH - mismatches - ignoredCharMatches)
				return true;
			else
			{
				ntMatches = 0; //must reset for each unsuccesful loop
				ignoredCharMatches = 0; //must reset for each unsuccesful loop
			}
		} //outer for() loop
		return false; // exhausted all possibilities, no matches found
	} //end isFuzzyMatch() method


	/**
	 * Combined isFuzzyMatch() and smithWaterman() to allow search that incorporates mismatches and/or indels
	 *
	 * @param query The query sequence to search for
	 * @param mismatches The maximum allowed number of mismatches between query and source
	 * @param ins The maximum allowed number of inserts in the SOURCE sequence (i.e. # deletions in query)
	 * @param del The maximum allowed number of deletions in the SOURCE sequence (i.e. # insertions in query)
	 * @return boolean TRUE if query is found within source, given the specified conditions; FALSE otherwise
	 * @see isFuzzyMatch()
	 * @see smithWaterman()
	 */
	public boolean isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int del)
	{
		Oligo swKey = new Oligo();

		try
		{
			swKey = smithWaterman(query, ins, del);
		}
		catch(OligoException error)
		{
			error.getStackTrace();
		}

		return isFuzzyMatch(swKey, mismatches);
	} //end isFuzzySmithWatermanMatch()


	public boolean isFuzzySmithWatermanMatch(Oligo query, int mismatches, int ins, int del, int minKeyLength)
	{
		Oligo swKey = new Oligo();

		try
		{
			swKey = smithWaterman(query, ins, del);
		}
		catch(OligoException error)
		{
			error.getStackTrace();
		}

		return isFuzzyMatch(swKey, mismatches) && swKey.length() >= minKeyLength;
	} //end isFuzzySmithWatermanMatch()


	/**
	 * Returns the last index of the first match, if found. Uses isFuzzyMatch() algorithm
	 *
	 * @param query The sequence to search for
	 * @param mismatches The max number of allowed mismatches
	 * @return int index (if found); otherwise, returns -1
	 */
	public int lastIndex(Oligo query, int mismatches)
	{
		return (isFuzzyMatch(query, mismatches)) ? getFirstMatchCoordinate(query, mismatches) + query.length() - 1 : -1;
	}


	/**
	 * Returns the length of oligo
	 *
	 * @return int - length of oligo
	 */
	public int length()
	{
		return oligo_length;
	}


	/**
	 * ligate the 5' (upstream) end of oligo2 to the 3' (downstream) end of oligo
	 *
	 * @param input
	 * @return
	 */
	public Oligo ligate(Oligo input)
	{
		return new Oligo(this.oligo + input.toString());
	}


	/**
	 * Induces base substitution mutations in oligo
	 *
	 * @param percent The percentage of mutation to be induced
	 * @return Oligo object - mutated oligo sequence
	 */
	public Oligo mutate(int percent)
	{
		int probability = (int) (Math.random() * (100 + 1)); //probability that the oligo will be mutated
		final int LAST_INDEX = oligo.length() - 1;
		StringBuffer tempOligo = new StringBuffer(oligo);

		if(probability > percent)
			return new Oligo(tempOligo);
		else if(probability <= percent)
		{
			int index = (int) (Math.random() * LAST_INDEX);
			char ch = generateRandomBase();
			if(oligo.charAt(index) != ch)
				tempOligo.setCharAt(index, ch);
		}
		return new Oligo(tempOligo);
	} //end mutate()


	/**
	 * Randomizes the given oligo sequence
	 *
	 * @return [Oligo object] randomized oligo sequence
	 */
	public Oligo randomize()
	{
		StringBuffer tempOligo = new StringBuffer(oligo);

		for(int i = 0; i <= 5 * oligo.length(); i++)
		{
			//randomly pick a character and move it to the end
			int randomIndex = (int) (Math.random() * (oligo.length() - 1));
			char tempChar = tempOligo.charAt(randomIndex);
			tempOligo.deleteCharAt(randomIndex);
			tempOligo.append(tempChar);
		}
		return new Oligo(tempOligo);
	} //end randomize()


	/**
	 * Resets the size of an oligo to zero length
	 *
	 * @return Oligo object: oligo of zero length
	 */
	public Oligo resetToNull()
	{
		return new Oligo();
	}


	/**
	 * Returns the reverse of the oligo
	 *
	 * @return Oligo object: the reversed oligo
	 */
	public Oligo reverse()
	{
		StringBuilder tempOligo = new StringBuilder(oligo);
		return new Oligo(tempOligo.reverse());
	}


	/**
	 * Sets the character to ignore during manipulation operations. Default is 'N'
	 *
	 * @param inputChar
	 */
	public void setIgnoredChar(char inputChar)
	{
		this.ignoredChar = inputChar;
	}


	/**
	 * Returns a modified query sequence based on success of a search. Allows user to specify maximum number of insertions and deletions in the source sequence
	 *
	 * @param query The query sequence to search for
	 * @param ins The maximum number of inserts allowed in the source sequence
	 * @param del The maximum number of deletes allowed in the source sequence
	 * @return Oligo object: transformed query sequence, if within specified conditions. Otherwise, returns "smithWaterman(): NotFoundException"
	 * @throws OligoException
	 */
	public Oligo smithWaterman(Oligo query, int ins, int del) throws OligoException
	{
		char[] s = ("x" + oligo).toUpperCase().toCharArray(); //pad 'x' as first char of oligo
		char[] q = ("x" + query).toUpperCase().toCharArray(); //pad 'x' as first char of query
		float[][] matrix = new float[q.length + 1][s.length + 1];
		float max = 0;
		float best = 0; //highest score in each comparison
		float iScore = 0;
		float jScore = 0;
		float diagScore = 0; //nucleotide-nucleotide alignment score: match = 1.0; mismatch = -0.3
		int imax = 0;
		int jmax = 0;

		//initialize row 0 to 0.0
		for(int i = 0; i <= q.length; i++)
			matrix[i][0] = (float) 0.0;

		//initialize column 0 to 0.0
		for(int j = 0; j <= s.length; j++)
			matrix[0][j] = (float) 0.0;

		//construct scores matrix
		for(int i = 1; i <= q.length - 1; i++)
		{
			for(int j = 1; j <= s.length - 1; j++)
			{
				//calculate diagonal score
				diagScore = (q[i] == s[j]) ? (float) 1.0 : (float) -0.3;
				best = matrix[i - 1][j - 1] + diagScore; //initially, asume diagonal is best score

				//calc max vertical score
				for(int vGap = i; vGap >= 0; vGap--)
				{
					iScore = matrix[i - vGap][j] - (float) (1 + 0.3 * vGap);
					best = (iScore > best) ? iScore : best;
				}

				//calc max horizontal score
				for(int hGap = j; hGap >= 0; hGap--)
				{
					jScore = matrix[i][j - hGap] - (float) (1 + 0.3 * hGap);
					best = (jScore > best) ? jScore : best;
				}

				matrix[i][j] = (best > 0) ? best : 0;
				if(matrix[i][j] > max)
				{
					max = matrix[i][j];
					imax = i;
					jmax = j;
				}
			} //for(int j...) loop
		} //for(int i...) loop

		//reconstruct alignment
		int nDel = 0;               //number of deletions in source
		int nIns = 0;               //number of insertions in source
		float above = 0;            //score of cell matrix[imax - 1][jmax]
		float left = 0;             //score of cell matrix[imax][jmax - 1]
		float diag = 0;             //score of cell matrix[imax - 1][jmax - 1]
		boolean isAbove = false;    //best score is cell above
		boolean isDiag = false;     //best score is diagonal cell
		boolean isLeft = false;     //best score is cell to left
		StringBuilder key = new StringBuilder();
		key = key.append(q[imax]);  //initialize key to ch associated with max score

		do
		{
			//System.out.println("imax = " + imax + " jmax = " + jmax + "  matrix[imax][jmax] = " + matrix[imax][jmax]);
			above = matrix[imax - 1][jmax];
			diag = matrix[imax - 1][jmax - 1];
			left = matrix[imax][jmax - 1];
			isAbove = (above >= diag) && (above >= left);
			isDiag = (diag >= above) && (diag >= left);
			isLeft = (left >= diag) && (left >= above);

			if(isAbove && !isDiag) //deletion in source sequence
			{
				imax--;
				nDel++;
			}
			else if(isDiag) //nucleotide-nucleotide match
			{
				imax--;
				jmax--;

				if(imax > 0) //prevent q[0] from being appended to key
					key.append(q[imax]);
			}
			else if(isLeft && !isDiag) //insertion in source sequence; insert 'n' into query to compensate
			{
				jmax--;
				key.append('n');
				nIns++;
			}
		} while(matrix[imax][jmax] > 0);

		//evaluate success of search based upon specified conditions
		if(nDel <= del && nIns <= ins)
			return new Oligo(key.reverse());
		else
			throw new OligoException(mnf, "smithWaterman()");
	} //end smithWaterman() method


	/**
	 * Deletes the first occurence of the target sequence based on a best fuzzy-match search. Overloaded to include the option of including indels
	 *
	 * @param query
	 * @param mismatches The maximum number of allowed mismatches
	 * @return Oligo object: spliced oligo, if target is found
	 * @throws OligoException
	 */
	public Oligo spliceOut(Oligo query, int mismatches) throws OligoException
	{
		if(isFuzzyMatch(query, mismatches))
		{
			int index = getFirstMatchCoordinate(query, mismatches);
			return new Oligo(oligo.substring(0, index) + oligo.substring(index + query.length(), oligo.length() - 1));
		}
		else
			throw new OligoException(mnf, "spliceOut()");
	} //end spliceOut()


	/**
	 * Splices out a given sequence based on specified start and end indices, inclusively
	 *
	 * @param start start index
	 * @param end end index
	 * @return Oligo object representing remaining sequence after splicing
	 * @throws OligoException
	 */
	public Oligo spliceOut(int start, int end) throws OligoException
	{
		String s = "";
		boolean startIndexOutofBounds = start < 0;
		boolean endIndexOutofBounds = end > oligo.length() - 1;

		if(start >= 0 && end <= oligo.length() - 1)
			s = oligo.substring(0, start) + oligo.substring(end + 1, oligo.length() - 1);
		else if(startIndexOutofBounds ^ endIndexOutofBounds)
		{
			if(startIndexOutofBounds)
				throw new OligoException(soob, "extractSequence()");
			else
				throw new OligoException(eoob, "extractSequence()");
		}
		else if(startIndexOutofBounds && endIndexOutofBounds)
			throw new OligoException("(Start&End)IndexOutofBounds", "extractSequence()");
		return new Oligo(s);
	} //end spliceOut()
} //end Oligo class
