// Fill out your copyright notice in the Description page of Project Settings.


#include "CombatSystem/Cell.h"
#include "CombatSystem/Figure.h"

// Sets default values
ACell::ACell()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	//PrimaryActorTick.bCanEverTick = true;

	// Dafault init
	State = ECellState::Empty;
	OccupiedFigure = nullptr;
	bIsHighlighted = false;
	Coordinates = FVector2D(0, 0);
	CellColor = FColor::White;

	// Create mesh-component (flat square)
	CellMesh = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("CellMesh"));
	RootComponent = CellMesh;
	// HERE set static mesh (ex. Plane) in constructor in UE Editor
}

// Called when the game starts or when spawned
void ACell::BeginPlay()
{
	Super::BeginPlay();
	// HERE additional initialization if needed
	
}

bool ACell::IsEmpty() const
{
	return State == ECellState::Empty && OccupiedFigure == nullptr;
}

void ACell::SetFigure(AFigure* newFigure)
{
	if (newFigure && State != ECellState::Inactive) {
		OccupiedFigure = newFigure;
		State = ECellState::Occupied;
		// Bijection: set a reference to a cell in a figure (NewFigure->SetCell(this); - add to Figure)
		newFigure->SetCell(this);
	}
}

void ACell::RemoveFigure()
{
	if (OccupiedFigure) { // != nullptr
		// Bijection: Clear the reference in the figure 
		OccupiedFigure->SetCell(nullptr);
		OccupiedFigure = nullptr;
	}
	State = ECellState::Empty;
}

void ACell::SetState(ECellState newState)
{
	State = newState;
	if (newState == ECellState::Empty) {
		RemoveFigure();
	}
}


// Called every frame
//void ACell::Tick(float DeltaTime)
//{
//	Super::Tick(DeltaTime);
//
//}

