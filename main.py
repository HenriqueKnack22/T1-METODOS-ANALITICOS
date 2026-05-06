"""
Simulador de Redes de Filas — Abordagem por Eventos Discretos

Implementa simulação de redes com topologia geral (incluindo ciclos e tandem),
roteamento estocástico entre filas e configuração variável de servidores.

Componentes principais (baseados no material da disciplina):
    - Fila:        estado e comportamento de uma fila G/G/c[/K]
    - Evento:      ocorrência com categoria, instante, fila de origem e destino
    - Escalonador: heap mínimo de eventos ordenados pelo instante de ocorrência

Geração de números pseudoaleatórios via LCG MINSTD / Park-Miller:
    X0 = semente, a = 16807, c = 0, M = 2^31 - 1
"""

import heapq
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import yaml

OUTPUT_ROOT = Path(__file__).parent / "saidas"
MODELOS_DIR = Path(__file__).parent / "cenarios"


# ===========================================================================
# Gerador LCG (MINSTD / Park-Miller)
# ===========================================================================

_estado_lcg = {"x": 7}

_A = 16807
_C = 0
_M = 2147483647


def proximo_aleatorio() -> float:
    """Avança o LCG e retorna um valor normalizado em (0, 1)."""
    _estado_lcg["x"] = (_A * _estado_lcg["x"] + _C) % _M
    return _estado_lcg["x"] / _M


def inicializar_semente(semente: int):
    _estado_lcg["x"] = semente


# ===========================================================================
# Fila
# ===========================================================================

class Fila:
    """
    Representa uma fila G/G/c[/K] com seus parâmetros e acumuladores.

    Atributos:
        servidores  — número de servidores paralelos
        capacidade  — limite máximo de clientes (None = sem limite)
        min_chegada — intervalo mínimo entre chegadas externas
        max_chegada — intervalo máximo entre chegadas externas
        min_servico — duração mínima de atendimento
        max_servico — duração máxima de atendimento
        na_fila     — clientes presentes no momento
        recusados   — clientes rejeitados por capacidade esgotada
        tempos      — tempo acumulado por estado (índice = número de clientes)
    """

    def __init__(
        self,
        nome: str,
        servidores: int,
        min_servico: float,
        max_servico: float,
        capacidade: int | None = None,
        min_chegada: float | None = None,
        max_chegada: float | None = None,
    ):
        self.nome = nome
        self.servidores = servidores
        self.capacidade = capacidade
        self.min_chegada = min_chegada
        self.max_chegada = max_chegada
        self.min_servico = min_servico
        self.max_servico = max_servico
        self.na_fila = 0
        self.recusados = 0
        if capacidade is None:
            self.tempos = [0.0]
        else:
            self.tempos = [0.0] * (capacidade + 1)

    def qtd(self) -> int:
        """Retorna o número atual de clientes na fila."""
        return self.na_fila

    def entrada(self):
        """Registra a entrada de um cliente."""
        self.na_fila += 1

    def saida_fila(self):
        """Registra a saída de um cliente."""
        self.na_fila -= 1

    def recusar(self):
        """Contabiliza um cliente recusado por fila cheia."""
        self.recusados += 1

    def aceita(self) -> bool:
        """True se há espaço para mais um cliente."""
        return self.capacidade is None or self.na_fila < self.capacidade

    def registrar_tempo(self, dt: float):
        """Acumula dt no estado atual; expande o vetor se necessário."""
        n = self.na_fila
        if n >= len(self.tempos):
            self.tempos.extend([0.0] * (n - len(self.tempos) + 1))
        self.tempos[n] += dt

    def notacao_kendall(self) -> str:
        """Retorna a notação G/G/c[/K] da fila."""
        if self.capacidade is None:
            return f"G/G/{self.servidores}"
        return f"G/G/{self.servidores}/{self.capacidade}"

    def distribuicao(self, tempo_total: float) -> list[float]:
        if tempo_total <= 0:
            return [0.0] * len(self.tempos)
        return [t / tempo_total for t in self.tempos]


# ===========================================================================
# Evento
# ===========================================================================

CHEGADA      = "CHEGADA"
PARTIDA      = "PARTIDA"
TRANSFERENCIA = "TRANSFERENCIA"


class Evento:
    """
    Representa um evento na simulação.

    categoria : CHEGADA | PARTIDA | TRANSFERENCIA
    instante  : tempo em que o evento ocorre
    origem    : fila de origem (PARTIDA e TRANSFERENCIA)
    destino   : fila de destino (CHEGADA e TRANSFERENCIA)
    """

    __slots__ = ("instante", "categoria", "origem", "destino", "_ordem")

    def __init__(
        self,
        instante: float,
        categoria: str,
        origem: str | None = None,
        destino: str | None = None,
    ):
        self.instante  = instante
        self.categoria = categoria
        self.origem    = origem
        self.destino   = destino
        self._ordem    = 0

    def __lt__(self, other: "Evento") -> bool:
        if self.instante != other.instante:
            return self.instante < other.instante
        return self._ordem < other._ordem


# ===========================================================================
# Escalonador
# ===========================================================================

class Escalonador:
    """Heap mínimo de eventos, ordenado pelo instante de ocorrência."""

    def __init__(self):
        self._eventos: list[Evento] = []
        self._contador = 0

    def inserir(self, evento: Evento):
        evento._ordem = self._contador
        self._contador += 1
        heapq.heappush(self._eventos, evento)

    def proximo_evento(self) -> Evento:
        return heapq.heappop(self._eventos)

    def esta_vazia(self) -> bool:
        return not self._eventos


# ===========================================================================
# Simulador
# ===========================================================================

def simular(
    filas: dict[str, Fila],
    rede: dict[str, list[tuple[float, str | None]]],
    chegadas_externas: dict[str, float],
    total_randoms: int,
    seed: int = 7,
) -> dict:
    """
    Executa a simulação da rede de filas pelo método de eventos discretos.

    Parâmetros:
        filas              — mapeamento nome → Fila
        rede               — mapeamento nome → [(probabilidade, destino|None)]
        chegadas_externas  — instante da primeira chegada por fila
        total_randoms      — quantidade de aleatórios a consumir (critério de parada)
        seed               — semente do LCG
    """
    inicializar_semente(seed)
    escalonador = Escalonador()

    relogio = 0.0
    restantes = [total_randoms]

    # --- Funções auxiliares ------------------------------------------------

    def intervalo_uniforme(minimo: float, maximo: float) -> float:
        r = proximo_aleatorio()
        restantes[0] -= 1
        return minimo + (maximo - minimo) * r

    def avanca_relogio(novo_tempo: float):
        nonlocal relogio
        dt = novo_tempo - relogio
        if dt > 0:
            for f in filas.values():
                f.registrar_tempo(dt)
        relogio = novo_tempo

    def escolher_destino(fila: Fila) -> str | None:
        """
        Sorteia o destino de um cliente ao terminar o serviço.
        Destino None = saída do sistema.
        Otimização: rota única determinística não consome número aleatório.
        """
        saidas = rede.get(fila.nome, [])
        if not saidas:
            return None
        if len(saidas) == 1 and abs(saidas[0][0] - 1.0) < 1e-12:
            return saidas[0][1]
        r = proximo_aleatorio()
        restantes[0] -= 1
        acumulado = 0.0
        for prob, destino in saidas:
            acumulado += prob
            if r < acumulado:
                return destino
        return None

    def planejar_servico(fila: Fila):
        destino = escolher_destino(fila)
        dt = intervalo_uniforme(fila.min_servico, fila.max_servico)
        if destino is None:
            ev = Evento(relogio + dt, PARTIDA, origem=fila.nome)
        else:
            ev = Evento(relogio + dt, TRANSFERENCIA, origem=fila.nome, destino=destino)
        escalonador.inserir(ev)

    def planejar_chegada(fila: Fila):
        dt = intervalo_uniforme(fila.min_chegada, fila.max_chegada)
        escalonador.inserir(Evento(relogio + dt, CHEGADA, destino=fila.nome))

    # --- Tratadores de eventos ---------------------------------------------

    def tratar_chegada(fila: Fila):
        if restantes[0] > 0:
            planejar_chegada(fila)
        if fila.aceita():
            fila.entrada()
            if fila.qtd() <= fila.servidores:
                if restantes[0] > 0:
                    planejar_servico(fila)
        else:
            fila.recusar()

    def tratar_partida(fila: Fila):
        fila.saida_fila()
        if fila.qtd() >= fila.servidores:
            if restantes[0] > 0:
                planejar_servico(fila)

    def tratar_transferencia(origem: Fila, destino: Fila):
        origem.saida_fila()
        if origem.qtd() >= origem.servidores:
            if restantes[0] > 0:
                planejar_servico(origem)
        if destino.aceita():
            destino.entrada()
            if destino.qtd() <= destino.servidores:
                if restantes[0] > 0:
                    planejar_servico(destino)
        else:
            destino.recusar()

    # --- Eventos iniciais --------------------------------------------------
    for nome, t_inicial in chegadas_externas.items():
        escalonador.inserir(Evento(t_inicial, CHEGADA, destino=nome))

    # --- Loop principal ----------------------------------------------------
    while restantes[0] > 0 and not escalonador.esta_vazia():
        ev = escalonador.proximo_evento()
        avanca_relogio(ev.instante)

        if ev.categoria == CHEGADA:
            tratar_chegada(filas[ev.destino])
        elif ev.categoria == PARTIDA:
            tratar_partida(filas[ev.origem])
        elif ev.categoria == TRANSFERENCIA:
            tratar_transferencia(filas[ev.origem], filas[ev.destino])

    # --- Drenagem de eventos remanescentes ---------------------------------
    while not escalonador.esta_vazia():
        ev = escalonador.proximo_evento()
        avanca_relogio(ev.instante)
        if ev.categoria == CHEGADA:
            f = filas[ev.destino]
            if f.aceita():
                f.entrada()
            else:
                f.recusar()
        elif ev.categoria == PARTIDA:
            filas[ev.origem].saida_fila()
        elif ev.categoria == TRANSFERENCIA:
            origem  = filas[ev.origem]
            destino = filas[ev.destino]
            origem.saida_fila()
            if destino.aceita():
                destino.entrada()
            else:
                destino.recusar()

    return {
        "tempo_global": relogio,
        "filas": filas,
    }


# ===========================================================================
# Saída de resultados
# ===========================================================================

def imprimir_resultados(resultado: dict, titulo: str):
    print("=" * 60)
    print(f"       SIMULADOR DE FILAS — {titulo}")
    print("=" * 60)
    print()
    print(f"Tempo total de simulação: {resultado['tempo_global']:.4f}")
    print()

    for fila in resultado["filas"].values():
        print(f"--- Fila {fila.nome}  ({fila.notacao_kendall()}) ---")
        if fila.min_chegada is not None:
            print(f"  Chegada:     {fila.min_chegada:.1f} .. {fila.max_chegada:.1f}")
        else:
            print(f"  Chegada:     (interna — alimentada por outras filas)")
        print(f"  Atendimento: {fila.min_servico:.1f} .. {fila.max_servico:.1f}")
        print(f"  Perdas:      {fila.recusados}")
        probs = fila.distribuicao(resultado["tempo_global"])
        print(f"  {'Estado':<8} {'Tempo':<18} {'Probabilidade':<15}")
        print("  " + "-" * 43)
        for i, (t, p) in enumerate(zip(fila.tempos, probs)):
            print(f"  {i:<8} {t:<18.4f} {p:<15.4f}")
        print("  " + "-" * 43)
        print(f"  {'Total':<8} {sum(fila.tempos):<18.4f} {sum(probs):<15.4f}")
        print()


def preparar_pasta(slug: str) -> Path:
    pasta = OUTPUT_ROOT / slug
    pasta.mkdir(parents=True, exist_ok=True)
    return pasta


def plotar_distribuicao(resultado: dict, titulo: str, arquivo: str | Path):
    """Gera gráfico de barras com a distribuição de probabilidade de cada fila."""
    filas = list(resultado["filas"].values())
    n = len(filas)
    fig, axes = plt.subplots(1, n, figsize=(max(6, 6 * n), 5), squeeze=False)

    for ax, fila in zip(axes[0], filas):
        probs  = fila.distribuicao(resultado["tempo_global"])
        estados = list(range(len(fila.tempos)))
        bars = ax.bar(estados, probs, color="steelblue", edgecolor="white", alpha=0.85)
        for bar, p in zip(bars, probs):
            if p > 0.005:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.005,
                    f"{p:.3f}",
                    ha="center", va="bottom", fontsize=8,
                )
        ax.set_xlabel("Estado (clientes)", fontsize=11)
        ax.set_ylabel("Probabilidade", fontsize=11)
        ax.set_title(f"Fila {fila.nome} — {fila.notacao_kendall()}", fontsize=12)
        if len(estados) <= 20:
            ax.set_xticks(estados)
        maxp = max(probs) if probs else 0
        ax.set_ylim(0, maxp * 1.15 if maxp > 0 else 1)
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle(titulo, fontsize=13)
    plt.tight_layout()
    plt.savefig(arquivo, dpi=150)
    plt.close(fig)
    print(f"Gráfico salvo em: {arquivo}")


def salvar_resultados(resultado: dict, titulo: str, arquivo: str | Path):
    with open(arquivo, "w", encoding="utf-8") as f:
        f.write(f"Simulador de Filas — {titulo}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Tempo total: {resultado['tempo_global']:.4f}\n\n")
        for fila in resultado["filas"].values():
            f.write(f"Fila {fila.nome}  ({fila.notacao_kendall()})\n")
            if fila.min_chegada is not None:
                f.write(f"  Chegada:     {fila.min_chegada:.1f} .. {fila.max_chegada:.1f}\n")
            else:
                f.write(f"  Chegada:     (interna)\n")
            f.write(f"  Atendimento: {fila.min_servico:.1f} .. {fila.max_servico:.1f}\n")
            f.write(f"  Perdas:      {fila.recusados}\n")
            probs = fila.distribuicao(resultado["tempo_global"])
            f.write(f"  {'Estado':<8} {'Tempo':<18} {'Probabilidade':<15}\n")
            f.write("  " + "-" * 43 + "\n")
            for i, (t, p) in enumerate(zip(fila.tempos, probs)):
                f.write(f"  {i:<8} {t:<18.4f} {p:<15.4f}\n")
            f.write("\n")
    print(f"Resultados salvos em: {arquivo}")


# ===========================================================================
# Carregamento de modelos YAML
# ===========================================================================
#
# Formato compatível com o simulador Java de referência:
#
#     rndnumbersPerSeed: 100000
#     seeds:
#     - 7
#
#     arrivals:
#        F1: 2.0           # instante da primeira chegada externa
#
#     queues:
#        F1:
#           servers: 1
#           capacity: 5    # omitir = capacidade infinita
#           minArrival: 2.0
#           maxArrival: 5.0
#           minService: 3.0
#           maxService: 5.0
#
#     network:
#     - source: F1
#       target: F2
#       probability: 0.8
#
# A probabilidade restante (1 − soma) é saída para o exterior.

def ler_modelo(yml_path: str | Path) -> dict:
    """
    Carrega um arquivo YAML e devolve as estruturas prontas para `simular`.
    """
    yml_path = Path(yml_path)
    texto = yml_path.read_text(encoding="utf-8")
    texto = "\n".join(linha for linha in texto.splitlines() if linha.strip() != "!PARAMETERS")
    dados = yaml.safe_load(texto) or {}

    filas: dict[str, Fila] = {}
    for nome, cfg in (dados.get("queues", {}) or {}).items():
        filas[nome] = Fila(
            nome=nome,
            servidores=int(cfg["servers"]),
            capacidade=cfg.get("capacity"),
            min_servico=float(cfg["minService"]),
            max_servico=float(cfg["maxService"]),
            min_chegada=(float(cfg["minArrival"]) if "minArrival" in cfg else None),
            max_chegada=(float(cfg["maxArrival"]) if "maxArrival" in cfg else None),
        )

    rede: dict[str, list[tuple[float, str | None]]] = {nome: [] for nome in filas}
    for entry in dados.get("network", []) or []:
        rede.setdefault(entry["source"], []).append((float(entry["probability"]), entry["target"]))

    chegadas_externas = {k: float(v) for k, v in (dados.get("arrivals", {}) or {}).items()}
    seeds = dados.get("seeds", [7]) or [7]

    return {
        "filas": filas,
        "rede": rede,
        "chegadas_externas": chegadas_externas,
        "total_randoms": int(dados.get("rndnumbersPerSeed", 100_000)),
        "seed": int(seeds[0]),
    }


# ===========================================================================
# Driver
# ===========================================================================

def rodar_simulacao(yml_path: Path):
    """Carrega o modelo, executa a simulação e grava os artefatos de saída."""
    nome = yml_path.stem

    print("=" * 60)
    print(f"  MODELO: {yml_path}")
    print("=" * 60)
    print()

    cfg = ler_modelo(yml_path)
    resultado = simular(
        filas=cfg["filas"],
        rede=cfg["rede"],
        chegadas_externas=cfg["chegadas_externas"],
        total_randoms=cfg["total_randoms"],
        seed=cfg["seed"],
    )

    imprimir_resultados(resultado, nome)

    pasta = preparar_pasta(nome)
    plotar_distribuicao(resultado, f"{nome} — Distribuição de Probabilidade", pasta / "distribuicao.png")
    salvar_resultados(resultado, nome, pasta / "resultados.txt")
    print()


def main():
    """
    Uso:
        python main.py                    # executa todos os .yml em cenarios/
        python main.py modelos/trab1.yml  # executa um modelo específico
        python main.py a.yml b.yml ...    # executa múltiplos modelos
    """
    args = sys.argv[1:]
    modelos = [Path(a) for a in args] if args else sorted(MODELOS_DIR.glob("*.yml"))

    if not modelos:
        print(f"Nenhum modelo .yml encontrado em {MODELOS_DIR}.")
        return

    for modelo in modelos:
        rodar_simulacao(modelo)


if __name__ == "__main__":
    main()
