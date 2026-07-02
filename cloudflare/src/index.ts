import { Container, getContainer } from "@cloudflare/containers";

// Streamlit is stateful (server-side session + WebSocket per browser tab),
// so every request is routed to a single default container instance rather
// than load-balanced across many — see docs/self-improvement-loop.md's
// sibling concern: this dashboard is a single-instance research tool, not a
// multi-tenant service. sleepAfter lets it scale to zero between visits.
export class StreamlitContainer extends Container {
  defaultPort = 8501;
  sleepAfter = "10m";

  override onStart() {
    console.log("Streamlit container started");
  }

  override onStop() {
    console.log("Streamlit container stopped");
  }

  override onError(error: unknown) {
    console.error("Streamlit container error:", error);
  }
}

interface Env {
  DASHBOARD_CONTAINER: DurableObjectNamespace<StreamlitContainer>;
}

export default {
  async fetch(request: Request, env: Env): Promise<Response> {
    // getContainer() with no name returns the same default instance for
    // every caller, which is what a stateful single-instance app needs.
    const container = getContainer(env.DASHBOARD_CONTAINER);
    return container.fetch(request);
  },
};
