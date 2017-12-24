
package jsat.linear;

/**
 * The value at a specified index for one dimension. This is a tool mean for use with sparce data structures. 
 * The values should not be backed by any list, and changes to the IndexValue should not alter any data 
 * structures. This class is mean to be returned by an iterator, and the iterator may reuse the same
 * IndexValue object for efficency. 
 * 
 * @author Edward Raff
 */
public final class IndexValue
{
    public int index;
    public double value;

    /**
     * Creates a new IndexValue
     * @param index the index for the given value
     * @param value the value at the specified index
     */
    public IndexValue(int index, double value)
    {
        this.index = index;
        this.value = value;
    }

}
