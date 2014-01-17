/*************************************************************************

 hv: main program

 ---------------------------------------------------------------------

                       Copyright (c) 2011
                Andreia P. Guerreiro <apg@dei.uc.pt>
                Carlos M. Fonseca <cmfonsec@dei.uc.pt>
                Michael T. M. Emmerich <emmerich@liacs.nl>

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA


*************************************************************************/

/*************************************************************************

 References:

 ---------------------------------------------------------------------
 [1] Andreia P. Guerreiro. Efficient algorithms for the assessment
    of stochastic multiobjective optimizers. Master’s thesis, IST,
    Technical University of Lisbon, Portugal, 2011.

 [2] Andreia P. Guerreiro, Carlos M. Fonseca, and Michael T. M. Emmerich.
    A fast dimension-sweep algorithm for the hypervolume indicator in four
    dimensions. In CCCG, pages 77–82, 2012.


*************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>


#if __GNUC__ >= 3
# define __hv_unused    __attribute__ ((unused))
#else
# define __hv_unused    /* no 'unused' attribute available */
#endif



typedef struct dlnode {
  double x[4];                    /* The data vector              */
  struct dlnode *nexty;
  struct dlnode *prevy;
  struct dlnode *nextz;
  struct dlnode *prevz;
  double xrightbelow; //x value of the closest point (according to x) with lower y and z values
} dlnode_t;


static int compare_node(const void *p1, const void* p2)
{
    const double x1 = *((*(const dlnode_t **)p1)->x+3);
    const double x2 = *((*(const dlnode_t **)p2)->x+3);

    return (x1 < x2) ? -1 : (x1 > x2) ? 1 : 0;
}

/*
 * Setup circular double-linked list in each dimension
 */


static dlnode_t *
setup_cdllist(double *data, int n, const double *ref)
{
    int d = 4;
    dlnode_t *head;
    dlnode_t **scratch;
    int i, j;
	double * data2;

    head  = malloc ((n+1) * sizeof(dlnode_t));



    for (i = 1; i <= n; i++) {
        for(j = 0; j < d; j++)
            head[i-1].x[j] = data[(i-1)*d + j];
        
    }
    
    head[n].x[0] = head[n].x[1] = head[n].x[2] = -DBL_MAX;
    head[n].x[3] = ref[3];


    scratch = malloc(n * sizeof(dlnode_t*));
        
    for (i = 0; i < n; i++){
        scratch[i] = head + i;
    }

    qsort(scratch, n, sizeof(dlnode_t*), compare_node);

    data2 = (double *) malloc (n * 4 * sizeof(double));
    for(i = 0; i < n; i++){
        for(j = 0; j < d; j++){
            data2[i*d + j] = scratch[i]->x[j];
        }
    }


    for(i = 0; i < n; i++){
        for(j = 0; j < d; j++)
            head[i].x[j] = data2[i*d + j];
    }

    free(data2);

    free(scratch);

    return head;
}

static void free_cdllist(dlnode_t * head)
{
    free(head);
}

static void free_sentinels(dlnode_t * list){

        free(list->nexty);
        free(list->prevy);

        free(list->nextz);
        free(list->prevz);
        
}

//box
typedef struct bx{
    double lx, ly, lz;
    double ux, uy, uz;

    struct bx * next;
    struct bx * prev;
}box;

//box list
typedef struct bxl{
    box * head;
    box * tail;
}boxlist;

//note: it is not necessary to define uz
static box * new_box(double lx, double ly, double lz, double ux, double uy, double uz){
    box * b = (box *)malloc(sizeof(box));
    
    b->lx = lx;
    b->ly = ly;
    b->lz = lz;
    
    b->ux = ux;
    b->uy = uy;
    b->uz = uz;
    
    return b;
}

//push_left == push_head
static void push_left(boxlist * bl, box * b){
    b->next = bl->head;
    bl->head->prev = b;
    b->prev = NULL;

    bl->head = b;
}


//may be empty
static void push_left_e(boxlist * bl, box * b){
    if(bl->head == NULL){
        bl->head = b;
        bl->tail = b;

        b->next = NULL;
        b->prev = NULL;

    }else{

        b->next = bl->head;
        bl->head->prev = b;
        b->prev = NULL;

        bl->head = b;
    }
}


//box b is freed and its volume is returned
static double close_box(box *b){
        double volume = (b->ux - b->lx) * (b->uy - b->ly) * (b->uz - b->lz);
        free(b);

        return volume;
}

static double shrink_box(box *b, double ux){
        double volume = (b->ux - ux) * (b->uy - b->ly) * (b->uz - b->lz);
        b->ux = ux;
        return volume;
}

static double close_boxes_right(boxlist *bl, double x, double uz){
    double volume = 0;	
    box * b = bl->tail;
    while(b->lx > x){
        bl->tail = b->prev;
        bl->tail->next = NULL;
        b->uz = uz;
        volume += close_box(b);
        b = bl->tail;
    }
    //b = bl->tail;
        if(x < b->ux){
        b->uz = uz;
        volume += shrink_box(b, x);
        
    }

    return volume;

}

static double close_boxes_left(boxlist *bl, double y, double uz){
    double volume;	
    box * b;

    double lastlx;
    double lx;
    double ly;
    double lz;

    if(y >= bl->head->uy)
        return 0;

    volume = 0;	
    b = bl->head;

    lastlx = b->ux;
    lx = b->lx;
    ly = b->ly;
    lz = b->lz;

    while(b != NULL && b->uy >= y){

        lastlx = b->ux;
        b->uz = uz;
        b->ly = y;
        bl->head = b->next;
        
        volume += close_box(b);
        b = bl->head;
    }

    if(bl->head != NULL)
        bl->head->prev = NULL;
    else
        bl->tail = NULL;

    b = new_box(lx, ly, lz, lastlx, y, 0);
    push_left_e(bl, b);

    return volume;

}


static double close_all_boxes(boxlist *bl, double uz){
    double volume = 0;
    box * b = bl->head;
    while(b != NULL){
        b->uz = uz;
        bl->head = b->next;
        volume += close_box(b);
        b = bl->head;
    }
    
    bl->head = bl->tail = NULL;
    return volume;
}



static void insert_after_y(dlnode_t *new, dlnode_t *prev){
    
    new->prevy = prev;
    new->nexty = prev->nexty;
    new->nexty->prevy = new;
    prev->nexty = new;

}

static void insert_after_z(dlnode_t *new, dlnode_t *prev){
    
    new->prevz = prev;
    new->nextz = prev->nextz;
    new->nextz->prevz = new;
    prev->nextz = new;

}


static void remove_y(dlnode_t *p){
    
    p->prevy->nexty = p->nexty;
    p->nexty->prevy = p->prevy;
    
}

static void remove_z(dlnode_t *p){
    
    p->prevz->nextz = p->nextz;
    p->nextz->prevz = p->prevz;
    
}


static boxlist * init_box_list(dlnode_t *p, dlnode_t * ynext){

    
    double previous_x = p->xrightbelow;

    dlnode_t * q = ynext; // y is swept in ascending order

    boxlist * bl = (boxlist *)malloc(sizeof(boxlist));
    box * b = (box *)malloc(sizeof(box));
    
    b->prev = b->next = NULL;
    bl->head = bl->tail = b; //box b (sentinel) will be removed later

    //box limits in x and y are defined (vertical boxes)
    while(q->x[0] > p->x[0] || q->x[2] > p->x[2]){

        // =
        if(q->x[2] <= p->x[2] && q->x[0] < previous_x && q->x[0] > p->x[0]){
            b = new_box(q->x[0], p->x[1], p->x[2], previous_x, q->x[1], 0);
            push_left(bl, b);

            previous_x = q->x[0];    
        }
        
        q = q->nexty;
    }
    
    b = new_box(p->x[0], p->x[1], p->x[2], previous_x, q->x[1], 0);
    push_left(bl, b);
    
    
    //the box created in the beginning is removed
    bl->tail = bl->tail->prev;
    free(bl->tail->next);
    bl->tail->next = NULL;

    
    return bl;
    
}


static double update_boxes(boxlist * bl, dlnode_t *p, dlnode_t * znext){
    
    double volume = 0;
    dlnode_t * q;
    
    q = znext;

    while(bl->head != NULL){      
            if(q->x[0] <= p->x[0]){

                if(q->x[1] <= p->x[1]){ // p is dominated in x and y

                    //close all boxes
                    //update volume                 
                    //last iteration
                     volume += close_all_boxes(bl, q->x[2]);
                }else{
                    
                    //update z
                    //close boxes in the left
                    //update volume
                     volume += close_boxes_left(bl, q->x[1], q->x[2]);
                }

                if(q->x[1] >= p->x[1] && p->x[0] < q->xrightbelow){
                    q->xrightbelow = p->x[0];
                }

                
            }else{ //since q does not dominate p and p does not dominate q, q has to
                    //to be to the right of p and below p
                

                //update z
                //closes boxes in the right
                //update volume
                 volume += close_boxes_right(bl, q->x[0], q->x[2]);

            }

        q = q->nextz;

    }

    free(bl);
    
    return volume;
    
}


//corresponds to HV4D - contribution function of the pseudocode
static double
hv_increment3DA(dlnode_t *p, dlnode_t * ynext, dlnode_t *znext){

     boxlist * bl = init_box_list(p, ynext);
//      if(bl == NULL)
//          return 0;

    return update_boxes(bl, p, znext);
}



//inner cycle (sweeps points according to coordinate z)
//note: returns the contribution of p only
//Function responsible for keeping points sorted according to coordinates y and z,
//for updating xrightbelow and for removing dominated points
static double
hv_increment3D(dlnode_t *list, dlnode_t *p, const double * ref){
    double loss = 0, gain;
    dlnode_t * q = list->prevz;
    dlnode_t * yprev = list->nexty;
    double xrightbelow = ref[0];
    dlnode_t * zprev =  list->nextz;
    dlnode_t * head3 = list->nextz;

    while(q != head3){

        if(q->x[2] >= p->x[2] && q->x[1] >= p->x[1] && q->x[0] >= p->x[0]){
            remove_y(q);
            remove_z(q);
            loss += hv_increment3DA(q, q->nexty, q->nextz);
            
        }else{

            if(q->x[2] > zprev->x[2] && (q->x[2] < p->x[2] || (q->x[2] == p->x[2] &&  q->x[1] <= p->x[1]))){
                 //because of xrightbelow, in order to avoid height 0 boxes (in spite of not causing wrong results)
                zprev = q;
            }

            if(q->x[1] > yprev->x[1] && (q->x[1] < p->x[1] || (q->x[1] == p->x[1] &&  q->x[2] <= p->x[2]))){ //testequal.4d.6
                yprev = q;
            }
            
            if(q->x[2] <= p->x[2] && q->x[1] <= p->x[1]){
 
              if(q->x[0] <= p->x[0]){
                return 0;
              }else if(q->x[0] < xrightbelow){ //q->x[0] > p->x[0]
                xrightbelow = q->x[0];
              }
            }
        }
        q = q->prevz;
    }


    p->xrightbelow = xrightbelow;


    gain = hv_increment3DA(p, yprev->nexty, zprev->nextz);

    //add p to the lists of coordinates two and three
    insert_after_z(p, zprev);
    insert_after_y(p, yprev);

    return gain - loss;
}


//outer cycle (sweeps list L_4)
static double
hv(dlnode_t *list, int c, const double * ref)
{

    dlnode_t *p = list;
    double hyperv = 0; //hvol4D
    double hypera = 0; //hvol3D
        
    int i;
    for (i = 0; i < c; i++) {

        hypera = hypera + hv_increment3D(list+c, p, ref);
        hyperv = hyperv + hypera * ((p+1)->x[3] - p->x[3]);
        p++;
    }

    return hyperv;
}






static void add_sentinels(dlnode_t * list, const double * ref){
    
    int i;
    
    for(i = 0; i <= 1; i++){
        dlnode_t * sentineli, * sentinelf;
        sentineli = (dlnode_t *) malloc(sizeof(dlnode_t));
        sentinelf = (dlnode_t *) malloc(sizeof(dlnode_t));

        sentineli->xrightbelow = -DBL_MAX;
        sentinelf->xrightbelow = -DBL_MAX;

        
        if(i == 0){
            sentineli->nexty = sentinelf;
            sentinelf->prevy = sentineli;
            list->nexty = sentineli;
            sentineli->prevy = list;
            list->prevy = sentinelf;
            sentinelf->nexty = list;
            
            
            sentineli->x[0] = ref[0];
            sentineli->x[1] = -DBL_MAX;
            sentineli->x[2] = -DBL_MAX;
            
            sentinelf->x[0] = -DBL_MAX;
            sentinelf->x[1] = ref[1];
            sentinelf->x[2] = -DBL_MAX;
        }else{
            sentineli->nextz = sentinelf;
            sentinelf->prevz = sentineli;
            list->nextz = sentineli;
            sentineli->prevz = list;
            list->prevz = sentinelf;
            sentinelf->nextz = list;
            
            sentineli->x[0] = DBL_MAX;
            sentineli->x[1] = DBL_MAX;
            sentineli->x[2] = -DBL_MAX;
            
            sentinelf->x[0] = -DBL_MAX;
            sentinelf->x[1] = -DBL_MAX;
            sentinelf->x[2] = ref[2];
            
        }

        sentineli->x[3] = DBL_MAX;
        sentinelf->x[3] = DBL_MAX;

        
    }
    
}


double guerreiro_hv4d(double *data, int n, const double *ref)
{
    dlnode_t *list;
    double hyperv;

    if (n == 0) { 
        /* Returning here would leak memory.  */
        hyperv = 0.0;
    } else {

        list = setup_cdllist(data, n, ref);
      
        add_sentinels(list+n, ref);

        hyperv = hv(list, n, ref);
        
        free_sentinels(list+n);
        /* Clean up.  */
        free_cdllist (list);
    
      
    }
    
    return hyperv;
}

