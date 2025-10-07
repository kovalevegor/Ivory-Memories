// Fill out your copyright notice in the Description page of Project Settings.


#include "CombatSystem/Figure.h"
#include "CombatSystem/Cell.h"

// Sets default values
AFigure::AFigure()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	//PrimaryActorTick.bCanEverTick = true;

	CurrentCell = nullptr;

}

void AFigure::SetCell(ACell* newCell)
{
	CurrentCell = newCell;
	if (newCell) {
		SetActorLocation(newCell->GetActorLocation() + FVector(0, 0, 10.0f)); // Elevate the figure above the cell
	}
}



// Called when the game starts or when spawned
//void AFigure::BeginPlay()
//{
//	Super::BeginPlay();
//	
//}
//
//// Called every frame
//void AFigure::Tick(float DeltaTime)
//{
//	Super::Tick(DeltaTime);
//
//}

