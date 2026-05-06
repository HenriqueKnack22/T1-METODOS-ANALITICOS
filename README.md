# Simulador de Redes de Filas — T1 Métodos Analíticos

Simulador de eventos discretos para redes de filas com notação G/G/c[/K]. Suporta topologias arbitrárias com roteamento probabilístico, realimentação (ciclos) e filas de capacidade infinita.

## Arquitetura

O simulador é composto por três entidades principais:

| Componente | Responsabilidade |
|---|---|
| `Fila` | Mantém o estado da fila (ocupação, perdas, tempos acumulados por estado) |
| `Evento` | Representa uma ocorrência tipada (`CHEGADA`, `PARTIDA`, `TRANSFERENCIA`) com instante e filas envolvidas |
| `Escalonador` | Heap mínimo que ordena eventos pelo instante de ocorrência |

A função `simular()` conduz o loop principal: extrai o próximo evento, avança o relógio global, acumula o tempo em cada estado e despacha para o tratador correspondente. O critério de parada é a quantidade de números aleatórios consumidos (`rndnumbersPerSeed`).

O gerador de números pseudoaleatórios é o LCG MINSTD / Park-Miller (`a = 16807`, `c = 0`, `M = 2³¹ − 1`).

## Configuração dos modelos (YAML)

Os cenários são descritos em arquivos `.yml` com formato compatível com o simulador Java de referência:

```yaml
rndnumbersPerSeed: 100000   # critério de parada

seeds:
- 7                          # semente do LCG (usa-se a primeira)

arrivals:
   F1: 2.0                   # instante da primeira chegada externa por fila

queues:
   F1:
      servers: 1
      capacity: 5             # omitir = capacidade infinita
      minArrival: 2.0
      maxArrival: 5.0
      minService: 3.0
      maxService: 5.0
   F2:
      servers: 2
      minService: 1.0
      maxService: 2.0

network:                      # roteamento entre filas
- source: F1
  target: F2
  probability: 0.8
- source: F2
  target: F1
  probability: 0.3
```

**Regras:**
- `capacity` omitido → fila com capacidade infinita (sem perdas por bloqueio).
- A probabilidade restante (`1 − soma das saídas`) representa saída para o exterior.
- Se `network` for omitido, todos os clientes saem do sistema após o atendimento.
- Ciclos e múltiplos destinos por fila são suportados.

## Como executar

Requer Python 3.12+ com `matplotlib` e `pyyaml` instalados.

```bash
python main.py                       # executa todos os .yml em cenarios/
python main.py cenarios/trab1.yml    # executa um modelo específico
python main.py a.yml b.yml ...       # executa múltiplos modelos
```

Os artefatos de saída são gravados em `saidas/<nome_do_modelo>/`:

```
saidas/
├── cenario1_gg1_5/     distribuicao.png + resultados.txt
├── cenario2_gg2_5/     distribuicao.png + resultados.txt
├── cenario3_tandem/    distribuicao.png + resultados.txt
├── cenario4_rede/      distribuicao.png + resultados.txt
└── trab1/              distribuicao.png + resultados.txt
```

## Cenários disponíveis

| Arquivo | Topologia | Descrição |
|---|---|---|
| `cenario1_gg1_5.yml` | G/G/1/5 | Fila simples — chegadas 2..5, atendimento 3..5 |
| `cenario2_gg2_5.yml` | G/G/2/5 | Dois servidores — chegadas 2..5, atendimento 3..5 |
| `cenario3_tandem.yml` | G/G/2/3 → G/G/1/5 | Rede em tandem com duas filas em série |
| `cenario4_rede.yml` | F1 → F2 / F3 / ext | Roteamento com três destinos possíveis |
| `trab1.yml` | F1 (∞) ↔ F2 (5) ↔ F3 (10) | Rede com realimentação entre três filas |

## Estrutura do repositório

```
.
├── main.py        # simulador completo (Fila, Evento, Escalonador, simular)
├── cenarios/      # cenários em YAML
├── referencia/    # simulador de referência do professor (validação)
└── saidas/        # saídas geradas em runtime
```
