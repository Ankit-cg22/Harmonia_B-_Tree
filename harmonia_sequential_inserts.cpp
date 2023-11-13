# include <bits/stdc++.h>
# define MAX 1000000000
using namespace std;

int k1=364; 
// k1= 3^0+...+3^5 = 364= sizeof(nodes)
// we assume 

int k2=121;   
//  k2=364-3^5= sizeof(prefix_sum) 
// we subtract size of leaf layer because leaf layer nodes do not have children and hence there wont be any entry for them in the prefix array

int m=2;            
// 1 key + 4 data items
 // so, each node can have max 3 keys and 4 children
int fanout=6;  
int nk=fanout-1;

bool flag=false;
int incremented_index=-1;

// B+ tree node 
typedef struct Node
{
  int* key;      // keys are sorted
  int** data;
}node;

node** nodes;
int* prefix_sum;

int lb(int* arr, int N, int X)
{
	int mid;

	// Initialise starting index and
	// ending index
	int low = 0;
	int high = N;

	// Till low is less than high
	while (low < high) {
		mid = low + (high - low) / 2;

		// If X is less than or equal
		// to arr[mid], then find in
		// left subarray
		if (X <= arr[mid]) {
			high = mid;
		}

		// If X is greater arr[mid]
		// then find in right subarray
		else {
			low = mid + 1;
		}
	}

	// Return the lb index
	return low;
}


void init(node* nn)
{

  // each node has nk keys , nk data arrays
  // each data arry has m integers 
  // here  nk= 5 , m = 2
  // means for each node there are 5 keys 
  // for each key we store an array of size 2 . the array contains 2 values : key , val
  nn->key=(int*)malloc(sizeof(int)*nk);                // init each node with 5 keys
  nn->data=(int**)malloc(sizeof(int*)*nk);             // and 5 data arrays, each having 2 integers
  for(int i=0;i<nk;i++)
  {
    nn->data[i]=(int*)malloc(sizeof(int)*m);
  }

  for(int i=0;i<nk;i++)
    nn->key[i]=MAX;

    for(int i=0;i<nk;i++)
    {
      for(int j=0;j<m;j++)
      {
        nn->data[i][j]=MAX;
      }
    }
}

int node_count=0; int prefix_sum_count=0;

int search(int key)
{
  int i=0; //node* ptr; ptr=nodes[0];
  int j=0;

  while(true)
  {
    //ptr=nodes[i];
    j=i;
    if(nodes[i]==NULL)                     // no node inserted yet , whole tree is empty
      return 0;
    int* key_arr=nodes[i]->key; // retrive the list of keys of this node 
    int ind=lb(key_arr,nk,key); // find lower bound of key in the array "key_arr" , ie find smallest value >= nk
    // we can do this because key_arr is sorted

    cout<<"i: "<<i<<"; ind: "<<ind<<" ; prefix_sum[i]: "<<prefix_sum[i]<<endl;

    if(prefix_sum[i]==MAX)                  // reached a leaf 
      break;

    i=prefix_sum[i]+ind; // using the formula child_idx = prefix_sum[node_idx] + ind (in research paper it was ind-1 , but here ind is already 0 based)

    if(key_arr[ind]==key)
      i++;

    // the right child of key contains the values >= key
    // so the entry for value = key will be present in right child(i currently pointing to left child , i+1 points to right child)
  }

  // j points to last non null node we reached

  int* key_arr=nodes[j]->key;
  int ind=lb(key_arr,nk,key); // we search for key in the key_array of nodes[j]
  //cout<<"SEARCH"<<endl;
  if(key_arr[ind]==key)
  {
    int** dd=nodes[j]->data;
    for(int r=0;r<m;r++)
        cout<<dd[ind][r]<<" "; // print the key value
    cout<<endl;
  }

  return j; // return the index where this key is present 
}

// print all the key value pairs where key is in range [k1 , k2]
void range_query(int k1,int k2)                   // scan leaves from first key to last key
{
  int j=search(k1); // find node where first key k1 is present
  int* key_arr=nodes[j]->key;
  int ind=lb(key_arr,nk,k1); // find the position of the key k1 in that node
  //cout<<"SEARCH"<<endl;
  cout<<"range_query: "<<endl;
  if(key_arr[ind]==k1)
  {
    // if k1 is present in this node , then print all the values from k1 till end(what ever is there , once we reach a key_arr[i] == MAX , it means that is end of keys of this node)
    cout<<"In"<<endl;
    ind++;
    int** dd=nodes[j]->data;

    for(;ind<nk;ind++)
    {
      if(key_arr[ind]==MAX)
        break;
      for(int r=0;r<m;r++)
          cout<<dd[ind][r]<<" ";
     }
    cout<<endl;
  }

  j++; // move to next node
  bool f2=false;

  while(nodes[j]!=NULL && f2==false)
  {
    int* key_arr=nodes[j]->key;
    int** dd=nodes[j]->data;
    int ind=0;
    for(;ind<nk ;ind++)
    {
      if(key_arr[ind]==MAX)
        break;
      if(dd[ind][0]>k2)
      {
        f2=true; // once we cross ahead of k2 we stop cause our range is [k1 , k2]
        break;
      }
      for(int r=0;r<m;r++)
          cout<<dd[ind][r]<<" ";
      cout<<endl;
    }
    j++;
  }
}

void put_in_middle(node* nn,int pos, int key, int* dd)   // put key and its dd at position pos in nn , nothing else . just put this key at position nn and move the keys from nn and ahead one step ahead
{
  // arr holds key of this node
  // and dt holds the data of this node
  int arr[nk]; int dt[nk][m];
  for(int i=0;i<nk;i++)
    arr[i]=nn->key[i];

  for(int i=0;i<nk;i++)
  {
    for(int j=0;j<m;j++)
    {
      dt[i][j]=nn->data[i][j];
    }
  }

  // we put the new key at its position
  // and move the keys ahead of it one step ahead
  nn->key[pos]=key;
  for(int i=pos+1;i<nk;i++)
    nn->key[i]=arr[i-1];

  for(int j=0;j<m;j++)
    nn->data[pos][j]=dd[j];
  // put the data of this node at pos 

  // move the data of the rest of the keys one step ahead
  for(int i=pos+1;i<nk;i++)
  {
    for(int j=0;j<m;j++)
    {
      nn->data[i][j]=dt[i-1][j];
    }
  }
}

void put_in_nodes(int ind,node* nn2)            // put nn2 in nodes at ind position
{
  node* g=nodes[ind]; nodes[ind]=nn2; int i=ind+1;
  // move all the nodes from ind+1 to k1 one position ahead 
  for(i=ind+1;i<k1 && nodes[i]!=NULL;i++)
  {
    node* h=nodes[i];
    nodes[i]=g;
    g=h;
  }
  nodes[i]=g;
}

void put_in_array(int* arr, int sz, int x, int ind)      // put value x at index ind in array arr. sz is size of array
{
  int g=arr[ind]; arr[ind]=x; int i=ind+1;
  // put x at position ind
  // move the old elements of index from ind to sz-1 one step ahead
  for(i=ind+1;i<k1 && arr[i]!=MAX;i++)
  {
    int h=arr[i];
    arr[i]=g;
    g=h;
  }
  arr[i]=g;
}

int parent(int ind)                           // give index of parent of node ind
{
  if(prefix_sum[0]==MAX || ind==0)                      // no element in prefix_sum
    return -1;

  int p=lb(prefix_sum,k2,ind); // find lower bound of ind in prefix_sum array (whose size = k2)
  //cout<<"pp: "<<p<<"; prefix_sum[p]: "<<prefix_sum[p]<<endl;
  if(prefix_sum[p]!=ind)
    p--;
  return p;
}


int insert_internal(node* nn,int ind,int key,int* dd,int orig)    // insert key in node* nn
{
  cout<<"start ind: "<<ind<<"; start key: "<<key<<endl;
  cout<<"nk: "<<nk<<endl;
  cout<<"Keys of start node nn: "<<endl;

  for(int i=0;i<nk;i++)
    cout<<nn->key[i]<<" ";
  cout<<endl;

  // find the position where we want this key to be inserted
  // key array is already sorted and we intend to keep it sorted
  // so we insert the new key such that key remains sorted
  int pos=lb(nn->key,nk,key);
  if(nn->key[nk-1]==MAX)                  // space there in node
  {
      cout<<"Inserting normally"<<endl;
      put_in_middle(nn,pos,key,dd); // this function simply puts this key at this specified position and move the keys present at pos and ahead one step ahead
      return 0;

  }
  else                                    // node is full
  {
    // node is full , need to create new nodes

    cout<<"Node is full"<<endl;
    for(int i=0;i<nk;i++)
    {
      cout<<nn->key[i]<<"# "<<endl;
    }
    cout<<"key: "<<key<<endl;

    int pos=lb(nn->key,nk,key); 
    node* nn2;  

    // middle_element : on basis of which we divide this current node into 2 parts
   int middle_element;  // middle_element goes to top next time
   int* dd_mid; 
   dd_mid=(int*)malloc(m*sizeof(int));

   cout<<"pos: "<<pos<<endl;

    if(prefix_sum[ind]==MAX)               // leaf
    {
      int mid=(nk+1)/2;
      cout<<"mid: "<<mid<<endl;

      if(pos<mid)
      {
        nn2=(node*)malloc(sizeof(node));
        init(nn2); 
        middle_element=nn->key[mid-1]; // mid-1 goes to top
        for(int j=0;j<m;j++)
          dd_mid[j]=nn->data[mid-1][j];

        // move all the keys from mid-1 to ahead to the new node nn2
        for(int i=mid-1;i<nk;i++)
        {
          nn2->key[i-(mid-1)]=nn->key[i];
        }

        for(int i=mid-1;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn2->data[i-(mid-1)][j]=nn->data[i][j];
          }
        }

        // now all the keys from mid-1 to nk are empty in the old node
        for(int i=mid-1;i<nk;i++)
        {
          nn->key[i]=MAX;
        }
        for(int i=mid-1;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }

        // now in the old node put the new key at the position we found 
        put_in_middle(nn,pos,key,dd);

      }
      else
      {
        // pos >=mid
        // we move the keys from mid to nk to new node
        // and the new key will be added to this new node
        if(pos==mid)
        {
          // the new key will be the key on basis of which we divide the 
          middle_element=key;
          for(int j=0;j<m;j++)
            dd_mid[j]=dd[j];
        }
        else
        {
          middle_element=nn->key[mid];
          cout<<"middle_element: "<<middle_element<<endl;
          for(int j=0;j<m;j++)
          {
            cout<<"j: "<<j<<endl;
            dd_mid[j]=nn->data[mid][j];
          }
        }

        // create and initialize new node
        nn2=(node*)malloc(sizeof(node));
        init(nn2);

        // move the keys from mid to nk to new node
        for(int i=mid;i<nk;i++)
        {
          nn2->key[i-mid]=nn->key[i];
        }
        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn2->data[i-mid][j]=nn->data[i][j];
          }
        }

        // put new key into the new node 
        put_in_middle(nn2,pos-mid,key,dd);


        // in the old node,  the keys from mid to nk are now empty , so set them to default ie empty slots
        // next when we want to enter new data , we can do in these empty slots
        cout<<"mid2: "<<mid<<endl;
        for(int i=mid;i<nk;i++)
        {
          cout<<"i: "<<i<<endl;
          nn->key[i]=MAX;
        }

        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }
      }

    }
    else // it is not a leaf
    {
      int mid=(nk+1)/2;

      if(pos<mid)
      {
        nn2=(node*)malloc(sizeof(node));
        init(nn2); 
        middle_element=nn->key[mid-1]; // mid-1 goes to top
        for(int j=0;j<m;j++)
          dd_mid[j]=nn->data[mid-1][j];

        // put mid to nk in new node
        for(int i=mid;i<nk;i++)                  // leaving mid-1
        {
          nn2->key[i-(mid)]=nn->key[i];
        }

        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn2->data[i-(mid)][j]=nn->data[i][j];
          }
        }

        // clear those positions from old node 
        for(int i=mid-1;i<nk;i++)
        {
          nn->key[i]=MAX;
        }
        for(int i=mid-1;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }

        // put nn in old node
        put_in_middle(nn,pos,key,dd);

      }
      else
      {
        if(pos==mid)
        {
          middle_element=key;
          for(int j=0;j<m;j++)
            dd_mid[j]=dd[j];
          nn2=(node*)malloc(sizeof(node));
          init(nn2);

          // if pos == mid we move mid to nk 
          // else we move from mid+1 to nk
          for(int i=mid;i<nk;i++)
          {
            nn2->key[i-mid]=nn->key[i];
          }
          for(int i=mid;i<nk;i++)
          {
            for(int j=0;j<m;j++)
            {
              nn2->data[i-mid][j]=nn->data[i][j];
            }
          }

        }
        else
        {
          // if pos > mid
          // then move mid+1 to nk to new node and add the new key to new node
          middle_element=nn->key[mid];
          for(int j=0;j<m;j++)
            dd_mid[j]=nn->data[mid][j];
          nn2=(node*)malloc(sizeof(node));
          init(nn2);

          for(int i=mid+1;i<nk;i++)
          {
            nn2->key[i-(mid+1)]=nn->key[i];
          }
          for(int i=mid+1;i<nk;i++)
          {
            for(int j=0;j<m;j++)
            {
              nn2->data[i-(mid+1)][j]=nn->data[i][j];
            }
          }
          put_in_middle(nn2,pos-(mid+1),key,dd);
        }

        for(int i=mid;i<nk;i++)
        {
          nn->key[i]=MAX;
        }
        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }
      }
    }

    // add the new node to the list of nodes
    put_in_nodes(ind+1,nn2);

    if(prefix_sum[ind]==MAX)               // leaf
    {
      int xx=0; int p=parent(ind);
      if(p==-1 || nodes[p]==NULL)     // leaf and parent is null means first node to be filled, and now, prefix_sum will get its first entry
      {
        cout<<endl;
        cout<<"Leaf with parent NULL:- ind: "<<ind<<"; p: "<<p<<endl;
        prefix_sum[0]=1;
        node* nn3; nn3=(node*)malloc(sizeof(node));
        init(nn3);
        //cout<<"middle_element: "<<middle_element<<endl;
        put_in_middle(nn3,0,middle_element,dd_mid); // put the mid element at index 0 of new node (the first one for this tree)
        //cout<<"middle_element: "<<middle_element<<endl;
        put_in_nodes(0,nn3); // put nn3 at index 0 as this is the first node in the tree
        return 0;
      }
      else
        xx=insert_internal(nodes[p],p,middle_element,dd_mid,ind);     // dd is NULL means we won't insert anything
      // recursively call this function to enter key = middle_element in the node stored at nodes[p] 
      // we move middle_element to parent , ie we move the value up into a new higher node 

      // xx means the index where the node was finally inserted
      ind+=xx;
      cout<<"Leaf with parent non-null:- ind+=xx: "<<ind<<"; p: "<<p<<endl;


      cout<<endl;
      p=parent(ind);
      int ind2=prefix_sum[p];
      cout<<"prefix_sum: "<<endl;
      for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
      {
        cout<<prefix_sum[i]<<" ";
      }
      cout<<endl;

      for(int i=0;i<k1;i++)
      {
        if(nodes[i]!=NULL)
        {
          cout<<"nodes["<<i<<"]->key: "<<endl;
          for(int j=0;j<nk;j++)
          {
            cout<<(nodes[i]->key)[j]<<" ";
          }
          cout<<endl;
          cout<<"nodes["<<i<<"]->data: "<<endl;
          for(int j=0;j<nk;j++)
          {
            cout<<"data["<<j<<"]- "<<endl;
            for(int f=0;f<m;f++)
            {
              cout<<(nodes[i]->data)[j][f]<<" ";
            }
            cout<<endl;
          }
          cout<<endl;
        }
      }

      cout<<"p: "<<p<<" ; ind: "<<ind<<"; ind2: "<<ind2<<endl;

      cout<<"prefix_sum in full leaf whose parent is not null, before: "<<endl;
      for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
      {
        cout<<prefix_sum[i]<<" ";
      }
      cout<<endl;

      //if(flag==false)
      //{
      cout<<"incremented_index: "<<incremented_index<<endl;
        for(int i=0;i<k2;i++)
        {
          if(prefix_sum[i]==MAX)
            break;

          if(prefix_sum[i]>ind2)
          {
            // all the nodes present ahead of prefix_sum[i] increase by 1 , cause the new index was inserted before them
            if(!(prefix_sum[i]==incremented_index && flag==true))
              prefix_sum[i]+=1;
          }
        }
      //}

      cout<<"prefix_sum in full leaf whose parent is not null,after : "<<endl;
      for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
      {
        cout<<prefix_sum[i]<<" ";
      }
      cout<<endl;

      return xx+1;
      // the current key was inserted at index xx+1 in nodes
    }
    else                                                    // non-leaf
    {
      int x=0; int no_of_keys_of_v1=0; int xx=0;
      for(int i=0;i<nk;i++)
      {
        if(nn->key[i]!=MAX)
          no_of_keys_of_v1++;
        else
            break;
      }
      //retriving index of parent node
      int p=parent(ind);
      //When parent node is null and insertion in non leaf node
      if(p==-1 || nodes[p]==NULL)                              // parent is NULL
      {
        cout<<endl;
        cout<<"In non-leaf with parent NULL, ind: "<<ind<<"; parent(ind): "<<p<<endl;
        //modification of prefix sum for the new parent node at 0th index
        put_in_array(prefix_sum,k2,1,0);
        x=2;
        //creation of new node (parent)
        node* nn3; nn3=(node*)malloc(sizeof(node));
        init(nn3);
        //cout<<"middle_element: "<<middle_element<<endl;
        //move the middle key of the child node in parent node
        put_in_middle(nn3,0,middle_element,dd_mid);
        //cout<<"middle_element: "<<middle_element<<endl;
        //put the newly created parent node in 0th Index of nodes array
        put_in_nodes(0,nn3);

        int pf_2=-1;
        int ss=prefix_sum[1]+2;
        cout<<"ss: "<<ss<<" ; orig+2: "<<orig+2<<"; ind: "<<ind<<endl;
        //division of key arrays based on middle element that has moved to new parent node
        //and based on number of keys, modifying the prefix sum array
        //calculating incremented index after insertion from root node
        if(ss>orig+2)
        {
          pf_2=ss-2+no_of_keys_of_v1;
          flag=false;
        }
        else
        {
          pf_2=ss-2+no_of_keys_of_v1+1;
          flag=true;
        }

        cout<<"pf_2: "<<pf_2<<endl;
        incremented_index=pf_2+2;;
        cout<<"incremented_index in orig+2: "<<incremented_index<<endl;
        //modyfying prefix sum for the second child of neww parent
        put_in_array(prefix_sum,k2,pf_2,2);

        cout<<"prefix_sum before "<<endl;
        for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
        {
          cout<<prefix_sum[i]<<" ";
        }
        cout<<endl;
       //modifying prefix sum 
        for(int i=1;i<k2;i++)
        {
          if(prefix_sum[i]==MAX)
            break;

          if(prefix_sum[i]>ind)
             // all the nodes present ahead of prefix_sum[i] increased by x(2) due to creation of 2 nodes (parent and 2nd child)
            prefix_sum[i]+=x;
        }
        
        cout<<"prefix_sum after: "<<endl;
        for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
        {
          cout<<prefix_sum[i]<<" ";
        }
        cout<<endl;
        //node was inserted at 2nd pos
        return 2;
      }
      else
      {
        //whwn parent is not null
        //find the index of parent node
        int p=parent(ind);
        int xx=insert_internal(nodes[p],p,middle_element,dd_mid,ind);
        //recursively call this function to enter key = middle_element in the node stored at nodes[p] 
        //we move middle_element to parent , ie we move the value up into a new higher node 
        ind+=xx;
        //xx - the index where the node was finally inserted
        x=1;
      }
      cout<<endl;
      cout<<"In non-leaf with parent non-null, ind+=xx: "<<ind<<"; parent(ind): "<<p<<endl;

      for(int i=0;i<k2;i++)
      {
        if(prefix_sum[i]==MAX)
          break;

        if(prefix_sum[i]>ind)
          // all the nodes present ahead of prefix_sum[i] increased by x(1) due to creation of node 
          prefix_sum[i]+=x;
      }

      //calculating incremented index after insertion 
      int pf_2=-1;
      int ss=prefix_sum[ind];
      cout<<"ss: "<<ss<<" ; orig+xx+1: "<<orig+xx+1<<"; ind: "<<ind<<"; xx: "<<xx<<"; orig: "<<orig<<endl;
      if(ss>orig+xx+1)
      {
        pf_2=ss+no_of_keys_of_v1;
        flag=false;
      }
      else
      {
        pf_2=ss+no_of_keys_of_v1+1;
        flag=true;
      }

      cout<<"pf_2: "<<pf_2<<endl;
      incremented_index=pf_2;
      cout<<"incremented_index in orig+ind: "<<incremented_index<<endl;
      //modyfing prefix sum for the another node created  after ind
      put_in_array(prefix_sum,k2,pf_2,ind+1);
      //the current key was inserted at index xx+1 in nodes
      return (xx+1);
    }

  }
}
                              
void insert(int* dd)
{
  // key = dd[0] , val = dd[1]
  int key=dd[0];
  int ind=search(key);             // gives index of leaf in nodes[], where key can be inserted
  node* nn=nodes[ind];

  if(nn==NULL)                      // nodes[0] is NULL
  {
    //new node is formed
    nn=(node*)malloc(sizeof(node));
    init(nn);
    nn->key[0]=key; // this is the only key in this node right now
    for(int j=0;j<m;j++)
    {
      nn->data[0][j]=dd[j];
    }

    nodes[0]=nn;
    // this is the only entry in this node right now , so we store it at 0th index
  }
  else                            // we have the leaf node nn to which original function can be applied
  {
    // we can insert this new key-val into one of the existing nodes
    insert_internal(nn,ind,key,dd,-1);
  }

  cout<<"prefix_sum: "<<endl;
  for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
  {
    cout<<prefix_sum[i]<<" ";
  }
  cout<<endl;

  // after inserting the key , for each key we print its kn keys and for each each key we print the 2 values : key , value
  cout<<"After inserting key "<<dd[0]<<": "<<endl;
  for(int i=0;i<k1;i++)
  {
    if(nodes[i]!=NULL)
    {
      cout<<"nodes["<<i<<"]->key: "<<endl;
      for(int j=0;j<nk;j++)
      {
        cout<<(nodes[i]->key)[j]<<" ";
      }
      cout<<endl;
      cout<<"nodes["<<i<<"]->data: "<<endl;
      for(int j=0;j<nk;j++)
      {
        cout<<"data["<<j<<"]- "<<endl;
        for(int f=0;f<m;f++)
        {
          cout<<(nodes[i]->data)[j][f]<<" ";
        }
        cout<<endl;
      }
      cout<<endl;
    }
  }

}

int main()
{
  // allocated memory for nodes(pointer to the tree)
  nodes=(node**)malloc(sizeof(node*)*k1);

  // allocate memory for the prefix sum tree
  prefix_sum=(int*)malloc(sizeof(int)*k2);

  // initialize all nodes as NULL nodes
  for(int i=0;i<k1;i++)nodes[i]=NULL;

  // initilize prefix sum with all values = MAX
  for(int i=0;i<k2;i++)prefix_sum[i]=MAX;
  
  int q;
  cout<<"Enter number of queries:- "<<endl;
  cin>> q;  // no_of_queries

  for(int i=0;i<q;i++)
  {
    cout<<"Enter 1 for inserting a tuple"<<endl;
    cout<<"Enter 2 for searching for a key"<<endl;
    cout<<"Enter 3 for doing range_query"<<endl;
    int option; cin>>option;
    if(option==1)                      // insert
    {
      int dd[m];
      cout<<"Enter "<<m<<" elements of a tuple, of which 1st should be key:- "<<endl;
      for(int i=0;i<m;i++)cin>>dd[i];
      flag=false;
      insert(dd);
    }
    else if(option==2)
    {
      cout<<"Enter key to be searched"<<endl; int ke; cin>>ke;
      search(ke);
    }
    else if(option==3)
    {
      cout<<"Enter 2 keys between which range_query is done"<<endl;
      int k1, k2; cin>>k1>>k2;
      range_query(k1,k2);
    }
  }

  return 0;
}
