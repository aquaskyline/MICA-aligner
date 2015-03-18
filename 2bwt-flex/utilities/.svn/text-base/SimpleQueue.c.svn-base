//
//    SimpleQueue.c
//
//    soap2-dp
//
//    Copyright (C) 2014, HKU
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 2
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
///////////////////////////////////////////////////////////////////////////////////////////////

#include "SimpleQueue.h"

SimpleQueue * SQCreate(int size) {
    SimpleQueue * ret = (SimpleQueue*) malloc (sizeof(SimpleQueue));
    ret->size = size;
    ret->count = 0;
    ret->head = 0;
    ret->tail = 0;
    ret->queueBody = (int*) malloc (sizeof(int) * ret->size);
    return ret;
}

void SQFree(SimpleQueue * queue) {
	if ( queue->queueBody!=NULL ) {
		free(queue->queueBody);
		queue->queueBody = NULL;
	}
    free(queue);
}

void SQInitialiseEmptyQueue(SimpleQueue * queue) {
    queue->count = 0;
    queue->head = 0;
    queue->tail = 0;
}

void SQInitialiseFullQueue(SimpleQueue * queue) {
    int i;
    SQInitialiseEmptyQueue(queue);
    for (i=0;i<queue->size;i++) {
        SQEnqueue(queue,i);
    }
}

int SQEnqueue(SimpleQueue * queue, int index) {
    if (queue->count == queue->size) {
        return 0;
    } else {
        queue->queueBody[queue->tail] = index;
        queue->tail = (queue->tail+1) % queue->size;
        queue->count++;
        return 1;
    }
}

int SQDequeue(SimpleQueue * queue, int * index) {
    if (queue->count == 0) {
        return 0;
    } else {
        (*index) = queue->queueBody[queue->head];
        queue->head = (queue->head+1) % queue->size;
        queue->count--;
        return 1;
    }
}

int SQIsFull(SimpleQueue * queue) {
    if (queue->count == queue->size) {
        return 1;
    }
    return 0;
}

int SQIsEmpty(SimpleQueue * queue) {
    if (queue->count == 0) {
        return 1;
    }
    return 0;
}

void SQPrintQueue(SimpleQueue * queue) {
    int i,j;
    printf("SQPrintQueue -- Print Queue Body\n");
    printf("SQInfo %d/%d used head@%d tail@%d\n",queue->count,queue->size,queue->head,queue->tail);
    i = 0;
    j = queue->head;
    for (;i<queue->count;i++) {
        printf("%5u*: %10u\n",i,queue->queueBody[j]);
        j = (j + 1) % queue->size;
    }
    for (;i<queue->size;i++) {
        printf("%5u : %10u\n",i,queue->queueBody[j]);
        j = (j + 1) % queue->size;
    }
}
