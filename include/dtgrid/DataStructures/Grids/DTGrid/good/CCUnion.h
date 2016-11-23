/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_ccunion_public_h
#define _grids_ccunion_public_h

namespace Grids
{

    /**
    * CCUnion_Public assumes that connected components are
    * pushed and popped like in a queue: first in first out
    */
    template<class Index, typename IndexArrayType = Index *> 
    class CCUnion_Public
    {
    public:
	struct CCI  // connected component interval
	{
	    unsigned int s;  // start
	    unsigned int e;  // end
	    unsigned int c;  // current
	    Index i;
	    bool start;      // start or end of connected component
	    CCI *next;
	};

	// Constants
	static const int CCI_MAX = 128;        // The maximal number of connected components allowed at a time
	static const int CCI_MASK = CCI_MAX-1;

    public:
	CCUnion_Public(Index G, IndexArrayType& indexPtr)
	    : indexPtr(indexPtr)
	{
	    this->G = G;
	    sp = 0;
	    ep = 0;
	}

	~CCUnion_Public()
	{
	}

	inline void addCCI(unsigned int s, unsigned int e)
	{
	    cci[ep].s = s;
	    cci[ep].e = e-1;
	    ep = (ep+1)&CCI_MASK;
	}

	inline void removeCCI()
	{
	    sp = (sp+1)&CCI_MASK;
	}

	inline void startCCUnion()
	{
	    unsigned int ccIterSp;
	    CCI *ptr;

	    if (sp==ep)
		return;

	    inside = 0;
	    done = false;
	    ccIterSp = sp;

	    cci[ccIterSp].c = cci[ccIterSp].s;
	    cci[ccIterSp].start = true;
	    cci[ccIterSp].i = getIndex(cci[ccIterSp].c) - G;
	    cci[ccIterSp].next = NULL;

	    ccs = &cci[ccIterSp];
	    ccIterSp = (ccIterSp+1)&CCI_MASK;

	    while(ccIterSp != ep)
	    {
		cci[ccIterSp].c = cci[ccIterSp].s;
		cci[ccIterSp].start = true;
		cci[ccIterSp].i = getIndex(cci[ccIterSp].c) - G;
		ptr = ccs;

		if (cci[ccIterSp].i < ptr->i)
		{
		    ccs = &cci[ccIterSp];
		    cci[ccIterSp].next = ptr;
		}
		else
		{
		    while(ptr->next != NULL && cci[ccIterSp].i>ptr->next->i)
		    {
			ptr = ptr->next;
		    }
		    // now cci[ccIterSp] should be inserted to the right of ptr
		    cci[ccIterSp].next = ptr->next;
		    ptr->next = &cci[ccIterSp];
		}

		ccIterSp = (ccIterSp+1)&CCI_MASK;
	    }
	}

	inline bool unionCCDone()
	{
	    return done;
	}


	inline void moveFirst()
	{
	    if (!( ccs->next==NULL || ccs->i < ccs->next->i ))
	    {
		CCI *mptr, *ptr;

		mptr = ccs;
		ptr = ccs = ccs->next;

		while(ptr->next != NULL && mptr->i > ptr->next->i)
		{
		    ptr = ptr->next;
		}
		// now cci[ccIterSp] should be inserted to the right of ptr
		mptr->next = ptr->next;
		ptr->next = mptr;
	    }
	}


	inline void getUnionCC(Index *i1, Index *i2)
	{

	    // WE KNOW THAT WE POINT TO THE START OF A CONNECTED COMPONENT
	    // AND THAT AT LEAST ONE CONNECTED COMPONENT EXISTS

	    inside = 1;
	    *i1 = ccs->i;
	    // find min
	    ccs->start = false;
	    ccs->c += 1;
	    ccs->i = getIndex(ccs->c) + G;
	    moveFirst();

	    while (inside > 0)
	    {
		// find min
		if (ccs->start)
		{
		    inside++;
		    ccs->start = false;
		    ccs->c += 1;
		    ccs->i = getIndex(ccs->c) + G;
		    moveFirst();
		}
		else
		{
		    inside--;
		    *i2 = ccs->i;
		    if (ccs->c == ccs->e)
		    {
			ccs = ccs->next;
		    }
		    else
		    {
			ccs->start = true;
			ccs->c += 1;
			ccs->i = getIndex(ccs->c) - G;
			moveFirst();
		    }


		    if (inside==0)
		    {
			// find min
			if (ccs==NULL)
			{
			    done = true;
			}
			else
			{
			    if (ccs->i <= *i2+1)
			    {
				inside++;
				ccs->start = false;
				ccs->c += 1;
				ccs->i = getIndex(ccs->c) + G;
				moveFirst();
			    }
			}
		    }

		}
	    }
	}

	Index getIndex(unsigned int i)
	{
	    return indexPtr[i];
	}

    protected:
	Index G;
	IndexArrayType& indexPtr;
	unsigned int inside;
	unsigned int sp, ep;   
	CCI cci[CCI_MAX];
	CCI *ccs;
	bool done;
	unsigned int numCCI;
    };


}
#endif
